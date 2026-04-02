/**
 * @file sledge_dev.c
 * @brief Search a protein sequence against a protein database using a custom pipeline.
 *
 * This program implements a custom version of HMMER's phmmer pipeline, optimized for
 * multi-threaded searches and sledgehammer-style filtering. It reads sequences in blocks,
 * assigns queries to worker threads, applies successive HMM filters and pipelines,
 * and stores or streams results based on user options.
 *
 */

#include <p7_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include <unistd.h>
#include "esl_threads.h"
#include "hmmer.h"
#include "sledge_dev.h"

void pipeline_thread(void *arg);

/**
 * @brief Print a standardized program banner.
 *
 * Outputs a header including the application name, SledgeHMMER version
 * and date, copyright, and license information.
 *
 * @param[in] fp        File stream to which the banner is printed (e.g., stdout).
 * @param[in] progname  Original program name (argv[0]) to extract the application name.
 * @param[in] banner    Short descriptive banner message.
 */
void
sledge_banner(FILE *fp, const char *progname, char *banner)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) esl_strdup(progname, -1, &appname);

  fprintf(fp, "# %s :: %s\n", appname, banner);
  fprintf(fp, "# SledgeHMMER %s (%s)\n", SHMMER_VERSION, SHMMER_DATE);
  fprintf(fp, "# %s\n", SHMMER_COPYRIGHT);
  fprintf(fp, "# %s\n", SHMMER_LICENSE);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname) free(appname);
  return;
}

/**
 * @brief Read a block of custom-formatted sequences into an ESL_SQ_BLOCK.
 *
 * This function reads up to sqBlock->listSize sequences from the input sequence file
 * into the provided ESL_SQ_BLOCK. It stops at EOF or when the block is full.
 *
 * @param[in]  sqfp      Sequence file pointer.
 * @param[out] sqBlock   Block structure to populate; its count and complete flag will be set.
 *
 * @return eslOK on success, eslEOF if no new sequences were read, eslENOSPACE if
 *         the block capacity is too small for all sequences in the file, eslEINVAL
 *         if listSize is not positive, or another Easel error code.
 */
 int
 read_cust(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock)
 {
   int status;
 
   sqBlock->count = 0;
 
   if (sqBlock->listSize <= 0)
     return eslEINVAL;
 
   while (sqBlock->count < sqBlock->listSize) {
     status = esl_sqio_Read(sqfp, sqBlock->list + sqBlock->count);
     if (status != eslOK)
       break;
     sqBlock->count++;
   }
 
   /* Treat EOF after at least one sequence as success for this read. */
   if (status == eslEOF && sqBlock->count > 0)
     status = eslOK;
 
   /* Block is full: if the sequence source is not exhausted, the block was too small. */
   if (status == eslOK && sqBlock->count == sqBlock->listSize) {
     if (sqfp->format == eslSQFILE_NCBI) {
       if ((uint32_t) sqfp->data.ncbi.index < sqfp->data.ncbi.num_seq)
         return eslENOSPACE;
     } else {
       FILE *fp = sqfp->data.ascii.fp;
       if (fp != NULL && !feof(fp))
         return eslENOSPACE;
     }
   }
 
   return status;
 }

/**
 * @brief Append a single sequence to an existing sequence block, growing it if needed.
 *
 * @param[in,out] db   Destination sequence block.
 * @param[in]     seq  Sequence to add.
 * @param[in]     abc  Alphabet object for block resizing.
 *
 * @return eslOK on success, or an error code on failure.
 */
int
add_seq(ESL_SQ_BLOCK *db, ESL_SQ *seq, const ESL_ALPHABET *abc)
{
  /* Check if db needs to be expanded */
  if (db->count == db->listSize)
    esl_sq_BlockGrowTo(db, db->listSize * 2, 1, abc);

  /* Add new sequence */
  esl_sq_Copy(seq, db->list + db->count);
  db->count++;

  return eslOK;
}

/**
 * @brief Append a new target hit to the result list, reallocating arrays if necessary.
 * 
 * Only appending failed hits. The default initialization for results is ACCEPT.
 *
 * @param[in,out] result     Result info struct to append to.
 * @param[in]     target_id  Identifier of the target sequence.
 * @param[in]     pid        Percent identity of the alignment.
 * @param[in]     eval       Log E-value of the hit.
 * @param[in]     length     Alignment length in residues.
 *
 * @return eslOK on success, eslENORESULT if memory allocation fails.
 */
static inline int
append_target(RESULT_INFO *result, const char *target_id, float pid, float eval, int length)
{
  /* Double capacity if needed */
  if (result->n_targets >= result->capacity) {
    int new_capacity = result->capacity * 2;
    char **new_targets = realloc(result->targets, new_capacity * sizeof(char*));
    float *new_pids    = realloc(result->pids,    new_capacity * sizeof(float));
    float *new_evals   = realloc(result->evals,   new_capacity * sizeof(float));
    int   *new_lengths = realloc(result->lengths, new_capacity * sizeof(int));
    if (!(new_targets && new_pids && new_evals && new_lengths)) {
      free(new_targets);
      free(new_pids);
      free(new_evals);
      free(new_lengths);
      return eslENORESULT;
    }
    result->targets  = new_targets;
    result->pids     = new_pids;
    result->evals    = new_evals;
    result->lengths  = new_lengths;
    result->capacity = new_capacity;
  }

  /* Add target information and mark as REJECT */
  result->targets[result->n_targets] = strdup(target_id);
  if (result->targets[result->n_targets] == NULL) {
    return eslENORESULT;
  }
  result->pids[result->n_targets]    = pid;
  result->evals[result->n_targets]   = eval;
  result->lengths[result->n_targets] = length;
  result->n_targets++;
  result->decision = REJECT;

  return eslOK;
}

/**
 * @brief Free all memory associated with a WORKER_INFO array across cores.
 *
 * @param[in,out] info       Array of WORKER_INFO to destroy.
 * @param[in]     num_cores  Number of worker threads.
 *
 * @return eslOK on success.
 */
int
destroy_info(WORKER_INFO *info, int num_cores)
{
  /* Loop through and free all members */
  for (int i = 0; i < num_cores; ++i) {
    for (int j = 0; j < info[i].num_queries; ++j) {
      for (int k = 0; k < info[i].result[j].n_targets; k++)
        free(info[i].result[j].targets[k]);
      free(info[i].result[j].targets);
      free(info[i].result[j].pids);
      free(info[i].result[j].evals);
      free(info[i].result[j].lengths);
      free(info[i].result[j].source);
    }
    free(info[i].result);
  }

  /* Now we can free the master object */
  free(info);

  return eslOK;
}

/**
 * @brief Buffer a line of output and flush when full.
 *
 * @param[in,out] buffer        Character buffer for output.
 * @param[in]     line          Line to append.
 * @param[in]     line_length   Length of the line string.
 * @param[in,out] si            SLEDGE_INFO containing output file and offset.
 */
void
add_to_buffer(char *buffer, char *line, size_t line_length, SLEDGE_INFO *si)
{
  /* If buffer will overflow, flush it */
  if (si->offset + line_length >= BUFFER_MAX) {
    fwrite_unlocked(buffer, 1, si->offset, si->out_fp);
    si->offset = 0;
  }

  /* Add to buffer */
  memcpy(buffer + si->offset, line, line_length);
  si->offset += line_length;
}

/**
 * @brief Aggregate and store results from worker threads into the output buffer.
 *
 * Depending on user --format mode, results presented accordingly.
 *
 * @param[in]  threadObj   Thread management object with per-thread data pointers.
 * @param[in]  si          Global SLEDGE_INFO with user options.
 * @param[out] results     Preallocated char buffer for result codes or lines.
 *
 * @return eslOK on success.
 */
int
store_results(ESL_THREADS *threadObj, SLEDGE_INFO *si, char *results)
{
  WORKER_INFO *info = NULL;
  int result_ctr = 0;
  char line[128];
  for (int i = 0; i < si->cores; i++) {
    info = (WORKER_INFO *)(threadObj->data[i]);
    for (int j = 0; j < info->num_queries; j++) {
      if (si->format == R_ALL) {
        /* Output ALL rejected hits for each query. */
        for (int hit = 0; hit < info->result[j].n_targets; hit++) {
          snprintf(line, sizeof(line), "%-15s %-15s %6.2f %7d %12.2E %8c\n",
                   info->result[j].source,
                   info->result[j].targets[hit],
                   info->result[j].pids[hit],
                   info->result[j].lengths[hit],
                   exp(info->result[j].evals[hit]),
                   info->result[j].decision);
          add_to_buffer(results, line, strlen(line), si);
        }
      }
      else if (si->format == R_REJECTQ) {
        /* Only print IDs of accepted sequences */
        if (info->result[j].decision == REJECT) {
          snprintf(line, sizeof(line), "%s\n", info->result[j].source);
          add_to_buffer(results, line, strlen(line), si);
        }
      }
      else if (si->format == R_REJECTT) {
        /* Only print IDs of failed comparisons */
        for (int hit = 0; hit < info->result[j].n_targets; hit++) {
          snprintf(line, sizeof(line), "%s\n", info->result[j].targets[hit]);
          add_to_buffer(results, line, strlen(line), si);
        }
      }
      else if (si->format == R_ACCEPT) {
        /* Only print IDs of accepted sequences */
        if (info->result[j].decision == ACCEPT) {
          snprintf(line, sizeof(line), "%s\n", info->result[j].source);
          add_to_buffer(results, line, strlen(line), si);
        }
      }
      else {
        /* Pack ACCEPT/REJECT decisions into a result string. */
        results[result_ctr++] = info->result[j].decision;
      }
    }
  }

  /* If it is an ACCEPT/REJECT string, terminate it */
  if (si->format == R_STRING)
    results[result_ctr] = '\0';

  return eslOK;
}

/**
 * @brief Get the current time in seconds using CLOCK_MONOTONIC_COARSE.
 *
 * Used for coarse timing of pipeline stages.
 *
 * @return Floating-point seconds since an arbitrary epoch.
 */
double
get_coarse_time()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC_COARSE, &ts);

  return ts.tv_sec + (ts.tv_nsec / 1.0e9);
}

/**
 * @brief Format a duration in seconds into an HH:MM:SS string.
 *
 * @param[in]  seconds       Time in seconds.
 * @param[out] buffer        Output buffer to receive formatted time.
 * @param[in]  buffer_size   Size of the output buffer.
 */
void
format_time(double seconds, char *buffer, size_t buffer_size)
{
  int hours   = (int)(seconds / 3600);
  int minutes = (int)((seconds - (hours * 3600)) / 60);
  int secs    = (int)(seconds - (hours * 3600) - (minutes * 60));
  snprintf(buffer, buffer_size, "%02d:%02d:%02d", hours, minutes, secs);
}

/**
 * @brief Print a dynamic progress bar with percentage, rate, elapsed, and ETA.
 *
 * @param[in] current Number of items processed.
 * @param[in] total   Total number of items.
 * @param[in] elapsed Elapsed time in seconds.
 */
void
print_progress(int current, int total, double elapsed)
{
  int   bar_width = 50;
  float progress  = (float)current / total;
  int   pos       = bar_width * progress;
  char  eta_str[10], total_str[10], bar[256];
  char *ptr = bar;

  /* Get ETA and elapsed time */
  format_time((elapsed / progress - elapsed), eta_str, sizeof(eta_str));
  format_time(elapsed, total_str, sizeof(total_str));

  /* Print progress bar in place */
  ptr += sprintf(ptr, "\rProgress: [");
  for (int i = 0; i < bar_width; i++)
    ptr += sprintf(ptr, "%s", i < pos ? "\xE2\x96\x88" : " ");

  sprintf(ptr, "] %.2f%% | %.2f it/s | Elapsed: %s | ETA: %s ",
          progress * 100,
          (elapsed == 0) ? 0.0 : (current / elapsed),
          total_str, eta_str);
  
  /* Add a new line if we reached the end */
  if (current == total) ptr += sprintf(ptr, "\n");
  fputs(bar, stderr);
}

/**
 * @brief Save a sequence block to a FASTA file with an indexed filename.
 *
 * Writes to a file named fpath_fid.fasta. If skip_lo < 0, writes every sequence;
 * otherwise skips indices in [skip_lo, skip_hi).
 *
 * @param[out] fp     File pointer (unused input; file is opened internally).
 * @param[in]  fpath  Base filepath prefix.
 * @param[in]  fid    File index to append to the prefix.
 * @param[in]  db     Sequence block to write.
 * @param[in]  skip_lo  Start of index range to omit, or -1 to write all.
 * @param[in]  skip_hi  End (exclusive) of index range to omit; ignored if skip_lo < 0.
 *
 * @return eslOK on success, or exits on error.
 */
int
save_seqs(FILE *fp, char *fpath, int fid, ESL_SQ_BLOCK *db, int skip_lo, int skip_hi)
{
  int    status = eslOK;
  char  *fl = NULL;
  int    len = snprintf(NULL, 0, "%s_%d.fasta", fpath, fid) + 1;

  ESL_ALLOC(fl, len * sizeof(char));
  snprintf(fl, len, "%s_%d.fasta", fpath, fid);

  if ((fp = fopen(fl, "w")) == NULL)
    p7_Fail("Failed to open file %s for writing\n", fl);

  for (int i = 0; i < db->count; i++) {
    if (skip_lo >= 0 && i >= skip_lo && i < skip_hi)
      continue;
    esl_sqio_Write(fp, db->list + i, eslSQFILE_FASTA, FALSE);
  }

  free(fl);
  fclose(fp);

  return status;

 ERROR:
  exit(status);
}

int
process_two_group_chunk(ESL_GETOPTS *go, ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *db_a, ESL_SQ_BLOCK *db_b,
                        ESL_SQ_BLOCK *check_a, ESL_SQ_BLOCK *check_b,
                        ESL_SQ_BLOCK *cdb, ESL_SQ_BLOCK *discard_db,
                        FILE *discard_fp, ESL_ALPHABET *abc, SLEDGE_INFO *si, bool assign_a,
                        char *results)
{
  int status = eslOK;
  int old_cores = si->cores;
  int old_qsize = si->qsize;
  ESL_SQ_BLOCK *chosen_db = assign_a ? db_a : db_b;
  ESL_SQ_BLOCK *other_db  = assign_a ? db_b : db_a;
  ESL_SQ_BLOCK *chosen_check = assign_a ? check_a : check_b;
  ESL_SQ_BLOCK *other_check  = assign_a ? check_b : check_a;

  if (chosen_check == NULL || chosen_check->count == 0) {
    memset(results, ACCEPT, qdb->count * sizeof(char));
    results[qdb->count] = '\0';
  } else {
    status = assign_master(go, qdb, chosen_check, si, results);
  }

  cdb->count = 0;
  discard_db->count = 0;
  for (int i = 0; i < qdb->count; i++) {
    if (results[i] == ACCEPT) add_seq(chosen_db, qdb->list + i, abc);
    else                      add_seq(cdb, qdb->list + i, abc);
  }

  if (cdb->count > 0) {
    if (si->cores > cdb->count) {
      si->cores = cdb->count;
      si->qsize = 1;
    } else {
      si->qsize = ceil((float) cdb->count / si->cores);
      si->cores = ceil((float) cdb->count / si->qsize);
    }

    if (other_check == NULL || other_check->count == 0) {
      memset(results, ACCEPT, cdb->count * sizeof(char));
      results[cdb->count] = '\0';
    } else {
      status = assign_master(go, cdb, other_check, si, results);
    }

    for (int i = 0; i < cdb->count; i++) {
      if (results[i] == ACCEPT) add_seq(other_db, cdb->list + i, abc);
      else                      add_seq(discard_db, cdb->list + i, abc);
    }
  }

  for (int i = 0; i < discard_db->count; i++)
    esl_sqio_Write(discard_fp, discard_db->list + i, eslSQFILE_FASTA, FALSE);

  si->cores = old_cores;
  si->qsize = old_qsize;
  return status;
}

/**
 * @brief Perform a three-way filtering workflow on a query block.
 *
 * Applies assign_master to db_2 and db_3, then partitions qdb into accepted,
 * rejected, and conditional blocks before reassigning deeper filters.
 *
 * @param[in] qdb         Query sequence block.
 * @param[in,out] db_1    Primary output block for triple-pass acceptances.
 * @param[in,out] db_2    Secondary block for mixed acceptances.
 * @param[in,out] db_3    Tertiary block for alternate acceptances.
 * @param[out] discard_fp File pointer to write discarded sequences.
 * @param[in,out] si      Global SLEDGE_INFO with user options.
 * @param[in] abc         Alphabet object for sequence blocks.
 * @param[in] go          Command line options.
 *
 * @return eslOK on success, or exits on memory or I/O error.
 */
int
check_triple(ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *db_1, ESL_SQ_BLOCK *db_2, ESL_SQ_BLOCK *db_3,
             FILE *discard_fp, SLEDGE_INFO *si, ESL_ALPHABET *abc, ESL_GETOPTS *go)
{
  ESL_SQ_BLOCK *cdb_2 = NULL;
  ESL_SQ_BLOCK *cdb_3 = NULL;
  char *results_2 = NULL;
  char *results_3 = NULL;
  int   status    = eslOK;
  int   tmp_ctr   = 0;

  /* Allocate result buffers sized to query count + 1 for terminator. */
  ESL_ALLOC(results_2, (qdb->count + 1) * sizeof(char));
  ESL_ALLOC(results_3, (qdb->count + 1) * sizeof(char));

  /* Get ACCEPT/REJECT decisions for db_1 */
  status = assign_master(go, qdb, db_2, si, results_2);
  status = assign_master(go, qdb, db_3, si, results_3);

  /* Create digital blocks for conditional candidates. */
  cdb_2 = esl_sq_CreateDigitalBlock(qdb->count, abc);
  cdb_3 = esl_sq_CreateDigitalBlock(qdb->count, abc);

  /* Save accept/candidate/discard sequences */
  for (tmp_ctr = 0; tmp_ctr < qdb->count; tmp_ctr++) {
    bool a2 = (results_2[tmp_ctr] == ACCEPT);
    bool a3 = (results_3[tmp_ctr] == ACCEPT);
    if (a2 && a3)                add_seq(db_1, qdb->list + tmp_ctr, abc);
    else if (!a2 && a3)          add_seq(cdb_2, qdb->list + tmp_ctr, abc);
    else if (a2 && !a3)          add_seq(cdb_3, qdb->list + tmp_ctr, abc);
    else                         esl_sqio_Write(discard_fp, qdb->list+tmp_ctr, eslSQFILE_FASTA, FALSE);
  }

  /* Store original cores and qsize before adjusting for smaller checks */
  int old_cores = si->cores, old_qsize = si->qsize;

  /* Process candidates in cdb_2 */
  if (cdb_2->count > 0) {
    /* Recompute cores and qsize */
    if (si->cores > cdb_2->count) { si->cores = cdb_2->count; si->qsize = 1; }
    else {
      si->qsize = ceil((float)cdb_2->count / si->cores);
      si->cores = ceil((float)cdb_2->count / si->qsize);
    }

    /* Get membership decision and write sequences accordingly */
    status = assign_master(go, cdb_2, db_1, si, results_2);
    for (int k = 0; k < cdb_2->count; k++) {
      if (results_2[k] == ACCEPT) add_seq(db_2, cdb_2->list+k, abc);
      else                        esl_sqio_Write(discard_fp, cdb_2->list+k, eslSQFILE_FASTA, FALSE);
    }
  }

  /* Restore original qsize and cores */
  si->cores = old_cores; si->qsize = old_qsize;

  /* Process candidates in cdb_3 */
  if (cdb_3->count > 0) {
    /* Recompute cores and qsize */
    if (si->cores > cdb_3->count) { si->cores = cdb_3->count; si->qsize = 1; }
    else {
      si->qsize = ceil((float)cdb_3->count / si->cores);
      si->cores = ceil((float)cdb_3->count / si->qsize);
    }

    /* Get membership decision and write sequences accordingly */
    status = assign_master(go, cdb_3, db_1, si, results_3);
    for (int k = 0; k < cdb_3->count; k++) {
      if (results_3[k] == ACCEPT) add_seq(db_3, cdb_3->list+k, abc);
      else                        esl_sqio_Write(discard_fp, cdb_3->list+k, eslSQFILE_FASTA, FALSE);
    }
  }

  /* Restore original qsize and cores */
  si->cores = old_cores; si->qsize = old_qsize;

  /* Free memory and destroy databases */
  esl_sq_DestroyBlock(cdb_2);
  esl_sq_DestroyBlock(cdb_3);
  free(results_2);
  free(results_3);

  return status;

 ERROR:
  exit(status);
}

/**
 * @brief Top-level assignment of queries to worker threads and collection of results.
 *
 * Creates WORKER_INFO structs, spawns threads, waits for completion, then stores results.
 *
 * @param[in] go    Command-line options.
 * @param[in] qdb   Query sequence block.
 * @param[in] sqdb  Target database sequence block.
 * @param[in,out] si   Global SLEDGE_INFO with options and settings.
 * @param[out] results Buffer to write per-query decision codes or E-values.
 *
 * @return eslOK on success, else an error code.
 */
int
assign_master(ESL_GETOPTS *go, ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *sqdb, SLEDGE_INFO *si, char *results)
{
  ESL_STOPWATCH *timer     = NULL;
  ESL_THREADS   *threadObj = NULL;
  WORKER_INFO   *info      = NULL;
  int            status    = eslOK;
  int            q_ctr     = 0;

  ESL_ALLOC(info, sizeof(*info) * si->cores);
  threadObj = esl_threads_Create(&pipeline_thread);
  
  /* Add db pointers and associated queries for each core */
  for (int i = 0; i < si->cores; i++) {
    info[i].si    = si;
    info[i].sqdb  = sqdb;
    info[i].qdb   = qdb;
    info[i].go    = go;
    info[i].start = q_ctr;
    if (si->qsize < qdb->count - q_ctr) {
      info[i].num_queries = si->qsize;
    } else {
      info[i].num_queries = qdb->count - q_ctr;
    }
    q_ctr += info[i].num_queries;
    ESL_ALLOC(info[i].result, info[i].num_queries * sizeof(RESULT_INFO));
    esl_threads_AddThread(threadObj, &info[i]);
  }

  /* Start and finish threads */
  esl_threads_WaitForStart(threadObj);
  esl_threads_WaitForFinish(threadObj);

  /* Process all results in requested format */
  status = store_results(threadObj, si, results);

  /* Free memory and destroy threads */
  destroy_info(info, si->cores);
  esl_threads_Destroy(threadObj);

  return status;

 ERROR:
  return status;
}

/**
 * @brief Worker thread entry point for processing query subsets.
 *
 * Implements early stopping where query is discarded after one failed comparison. 
 * Each thread builds HMM models for its assigned queries and scans the target database,
 * applying custom pipeline logic (p7_Pipeline_cust) for each sequence.
 *
 * @param[in,out] arg  Pointer to ESL_THREADS object containing per-thread data.
 */
void
pipeline_thread(void *arg)
{
  int            i, q_ctr, status, workeridx, seed;
  ESL_ALPHABET  *abc;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  P7_BUILDER    *bld = NULL;
  P7_BG         *bg = NULL;

  /* Initialize alphabet and start this thread */
  impl_Init();
  abc = esl_alphabet_Create(eslAMINO);
  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  /* Fetch info object for this thread */
  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);
  ESL_GETOPTS *go = info->go;
  bg = p7_bg_Create(abc);

  /* Initialize p7 builder */
  bld = p7_builder_Create(NULL, abc);
  if ((seed = esl_opt_GetInteger(go, "--seed")) > 0) {
    esl_randomness_Init(bld->r, seed);
    bld->do_reseeding = TRUE;
  }
  bld->EmL = esl_opt_GetInteger(go, "--EmL");
  bld->EmN = esl_opt_GetInteger(go, "--EmN");
  bld->EvL = esl_opt_GetInteger(go, "--EvL");
  bld->EvN = esl_opt_GetInteger(go, "--EvN");
  bld->EfL = esl_opt_GetInteger(go, "--EfL");
  bld->EfN = esl_opt_GetInteger(go, "--EfN");
  bld->Eft = esl_opt_GetReal(go, "--Eft");

  /* Load score system */
  status = p7_builder_LoadScoreSystem(bld,
            esl_opt_GetString(go, "--mx"),
            esl_opt_GetReal(go, "--popen"),
            esl_opt_GetReal(go, "--pextend"),
            bg);
  if (status != eslOK)
    p7_Fail("Failed to set single query seq score system:\n%s\n", bld->errbuf);

  /* Check each query sequence for similarity against target db */
  for (q_ctr = 0; q_ctr < info->num_queries; q_ctr++) {
    P7_OPROFILE *om    = NULL;
    P7_TOPHITS  *th    = p7_tophits_Create();
    P7_PIPELINE *pli   = NULL;
    ESL_SQ      *qsq   = info->qdb->list + info->start + q_ctr;

    /* Initialize result structure */
    RESULT_INFO *ri = &info->result[q_ctr];
    ri->n_targets = 0;
    ri->max_pid   = -1.0;
    ri->decision  = ACCEPT;
    ri->capacity  = INITIAL_HITS_CAPACITY;
    ri->source    = strdup(qsq->name);
    ri->targets   = malloc(INITIAL_HITS_CAPACITY * sizeof(char*));
    ri->pids      = malloc(INITIAL_HITS_CAPACITY * sizeof(float));
    ri->evals     = malloc(INITIAL_HITS_CAPACITY * sizeof(float));
    ri->lengths   = malloc(INITIAL_HITS_CAPACITY * sizeof(int));

    /* Build query profile */
    p7_SingleBuilder(bld, qsq, bg, NULL, NULL, NULL, &om);
    pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS);
    p7_pli_NewModel(pli, om, bg);

    /* Check against each target sequence */
    for (i = 0; i < info->sqdb->count; i++) {
      ESL_SQ *tsq = info->sqdb->list + i;
      p7_pli_NewSeq(pli, tsq);
      p7_bg_SetLength(bg, tsq->n);
      p7_oprofile_ReconfigLength(om, tsq->n);

      /* Run the p7 pipeline against current pair of sequences */
      status = p7_Pipeline_cust(pli, om, bg, tsq, qsq, th, ri, info->si);
      if (ri->decision == REJECT && info->si->early_stop) break;
    }

    /* Destroy p7 objects */
    p7_tophits_Destroy(th);
    p7_pipeline_Destroy(pli);
    p7_oprofile_Destroy(om);
  }

  /* Finish this thread and free remaining p7 objects */
  esl_threads_Finished(obj, workeridx);
  p7_bg_Destroy(bg);
  p7_builder_Destroy(bld);
}

/**
 * @brief Custom pipeline that runs heuristic filters, HMM alignment to compute similarity score.
 * 
 * @param[in]  pli     Pipeline object with filters and thresholds.
 * @param[in]  om      Profile HMM object for the query.
 * @param[in]  bg      Background model object.
 * @param[in]  sq      Target sequence to scan.
 * @param[in]  qsq     Query sequence (for PID calculation).
 * @param[in]  hitlist Top hits structure (unused here).
 * @param[in]  result  Result structure to append hits and decision.
 * @param[in]  si      Global SLEDGE_INFO with thresholds and options.
 *
 * @return eslOK on success, or an error code on failure.
 */
int
p7_Pipeline_cust(P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg,
                 const ESL_SQ *sq, const ESL_SQ *qsq,
                 P7_TOPHITS *hitlist, RESULT_INFO *result,
                 SLEDGE_INFO *si)
{
  float usc, vfsc, fwdsc, filtersc, nullsc, seq_score;
  double P;
  int status;

  if (sq->n == 0) return eslOK;
  if (sq->n > 100000)
    ESL_EXCEPTION(eslETYPE,
      "Target sequence length >100K, over pipeline limit.");

  p7_omx_GrowTo(pli->oxf, om->M, 0, sq->n);
  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

  /* MSV filter */
  p7_MSVFilter(sq->dsq, sq->n, om, pli->oxf, &usc);
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) {
    if (result->n_targets == 0) result->decision = ACCEPT;
    return eslOK;
  }
  pli->n_past_msv++;

  /* Bias filter */
  if (pli->do_biasfilter) {
    p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
    seq_score = (usc - filtersc) / eslCONST_LOG2;
    P = esl_gumbel_surv(seq_score, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
    if (P > pli->F1) {
      if (result->n_targets == 0) result->decision = ACCEPT;
      return eslOK;
    }
  } else filtersc = nullsc;
  pli->n_past_bias++;

  /* Scan-specific setup */
  if (pli->mode == p7_SCAN_MODELS) {
    if (pli->hfp) p7_oprofile_ReadRest(pli->hfp, om);
    p7_oprofile_ReconfigRestLength(om, sq->n);
    if ((status = p7_pli_NewModelThresholds(pli, om)) != eslOK)
      return status;
  }

  /* Viterbi filter */
  if (P > pli->F2) {
    p7_ViterbiFilter(sq->dsq, sq->n, om, pli->oxf, &vfsc);
    seq_score = (vfsc - filtersc) / eslCONST_LOG2;
    P = esl_gumbel_surv(seq_score, om->evparam[p7_VMU], om->evparam[p7_VLAMBDA]);
    if (P > pli->F2) {
      if (result->n_targets == 0) result->decision = ACCEPT;
      return eslOK;
    }
  }
  pli->n_past_vit++;

  /* Forward parser */
  p7_ForwardParser(sq->dsq, sq->n, om, pli->oxf, &fwdsc);
  seq_score = (fwdsc - filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  if (P > pli->F3) {
    if (result->n_targets == 0) result->decision = ACCEPT;
    return eslOK;
  }
  pli->n_past_fwd++;

  /* Backward and domain definition */
  p7_omx_GrowTo(pli->oxb, om->M, 0, sq->n);
  p7_BackwardParser(sq->dsq, sq->n, om, pli->oxf, pli->oxb, NULL);
  return p7_domaindef_cust(sq, qsq, om, pli->oxf, pli->oxb,
                           pli->fwd, pli->bck, pli->ddef, bg,
                           FALSE, NULL, NULL, NULL,
                           result, &nullsc, si);
}

/**
 * @brief Domain definition and clustering custom workflow.
 *
 * Performs posterior decoding, null2 recalculation, and cluster-based envelope resolution.
 *
 * @param[in]  sq      Target sequence.
 * @param[in]  qsq     Query sequence.
 * @param[in]  om      Profile HMM.
 * @param[in]  oxf     Forward scores matrix.
 * @param[in]  oxb     Backward scores matrix.
 * @param[in]  fwd     Forward matrix for null2.
 * @param[in]  bck     Backtrack matrix for null2.
 * @param[in]  ddef    Domain definition object.
 * @param[in]  bg      Background model.
 * @param[in]  long_target  Flag for long target pipelines.
 * @param[in]  bg_tmp  Temporary background for null2.
 * @param[in]  scores_arr   Scratch array for domain scoring.
 * @param[in]  fwd_emissions_arr Scratch array for emissions.
 * @param[in]  result  Result struct to append domains.
 * @param[in]  nullsc  Null model score pointer.
 * @param[in]  si      Global SLEDGE_INFO.
 *
 * @return eslOK on success, or eslENORESULT on growth or decoding error.
 */
int
p7_domaindef_cust(const ESL_SQ *sq, const ESL_SQ *qsq, P7_OPROFILE *om,
                  P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck,
                  P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
                  P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr,
                  RESULT_INFO *result, float *nullsc, SLEDGE_INFO *si)
{
   int i, j;
   int triggered;
   int d;
   int i2,j2;
   int last_j2;
   int nc;
   int saveL     = om->L;	                          /* Save the length config of <om>; will restore upon return                               */
   int save_mode = om->mode;	                        /* Likewise for the mode.                                                                 */
   int status;
   
   if ((status = p7_domaindef_GrowTo(ddef, sq->n))      != eslOK) return eslENORESULT;  /* ddef's btot,etot,mocc now ready for seq of length n */
   if ((status = p7_DomainDecoding(om, oxf, oxb, ddef)) != eslOK) return eslENORESULT;  /* ddef->{btot,etot,mocc} now made.                    */
 
   esl_vec_FSet(ddef->n2sc, sq->n+1, 0.0);          /* ddef->n2sc null2 scores are initialized                                                 */
   ddef->nexpected = ddef->btot[sq->n];             /* posterior expectation for # of domains (same as etot[sq->n])                            */
 
   p7_oprofile_ReconfigUnihit(om, saveL);	          /* process each domain in unihit mode, regardless of om->mode                              */
   i = -1;
   triggered = FALSE;
 
   for (j = 1; j <= sq->n; j++) {
 
     if (! triggered) {
       /* xref J2/101 for what the logic below is: */
       if       (ddef->mocc[j] - (ddef->btot[j] - ddef->btot[j-1]) <  ddef->rt2) i = j;
       else if  (i == -1)                                                        i = j;
       if       (ddef->mocc[j]                                     >= ddef->rt1) triggered = TRUE;
     }
     else if (ddef->mocc[j] - (ddef->etot[j] - ddef->etot[j-1])  <  ddef->rt2) {
         /* We have a region i..j to evaluate. */
         p7_omx_GrowTo(fwd, om->M, j-i+1, j-i+1);
         p7_omx_GrowTo(bck, om->M, j-i+1, j-i+1);
         ddef->nregions++;
         if (is_multidomain_region(ddef, i, j)) {
             /* This region appears to contain more than one domain, so we have to
              * resolve it by cluster analysis of posterior trace samples, to define
              * one or more domain envelopes.
              */
             ddef->nclustered++;
 
             /* Resolve the region into domains by stochastic trace
              * clustering; assign position-specific null2 model by
              * stochastic trace clustering; there is redundancy
              * here; we will consolidate later if null2 strategy
              * works
              */
             p7_oprofile_ReconfigMultihit(om, saveL);
             p7_Forward(sq->dsq+i-1, j-i+1, om, fwd, NULL);
 
             region_trace_ensemble(ddef, om, sq->dsq, i, j, fwd, bck, &nc);
             p7_oprofile_ReconfigUnihit(om, saveL);
             /* ddef->n2sc is now set on i..j by the traceback-dependent method */
 
             last_j2 = 0;
             for (d = 0; d < nc; d++) {
                   p7_spensemble_GetClusterCoords(ddef->sp, d, &i2, &j2, NULL, NULL, NULL);
                   if (i2 <= last_j2) ddef->noverlaps++;
 
                   /* Note that k..m coords on model are available, but
                      * we're currently ignoring them.  This leads to a
                      * rare clustering bug that we eventually need to fix
                      * properly [xref J3/32]: two different regions in one
                      * profile HMM might have hit same seq domain, and
                      * when we now go to calculate an OA trace, nothing
                      * constrains us to find the two different alignments
                      * to the HMM; in fact, because OA is optimal, we'll
                      * find one and the *same* alignment, leading to an
                      * apparent duplicate alignment in the output.
                      *
                      * Registered as #h74, Dec 2009, after EBI finds and
                      * reports it.  #h74 is worked around in p7_tophits.c
                      * by hiding all but one envelope with an identical
                      * alignment, in the rare event that this
                      * happens. [xref J5/130].
                   */
                   ddef->nenvelopes++;
 
                   /*the !long_target argument will cause the function to recompute null2
                    * scores if this is part of a long_target (nhmmer) pipeline */
                   if ((status = rescore_domain_cust(ddef, om, sq, qsq, fwd, bck, i2, j2, TRUE, bg, long_target, bg_tmp, scores_arr, fwd_emissions_arr, result, nullsc, si)) == eslOK)
                     last_j2 = j2;
                   if ((result->decision == REJECT) && (si->early_stop)) break; // breaking out of first loop over domains if early stopping
 
             }
 
             p7_spensemble_Reuse(ddef->sp);
             p7_trace_Reuse(ddef->tr);
             if ((result->decision == REJECT) && (si->early_stop)) break; // must break out of second loop if early stopping!!
         }
         else {
             /* The region looks simple, single domain; convert the region to an envelope. */
             ddef->nenvelopes++;
             status = rescore_domain_cust(ddef, om, sq, qsq, fwd, bck, i, j, FALSE, bg, long_target, bg_tmp, scores_arr, fwd_emissions_arr, result, nullsc, si);
             if ((result->decision == REJECT) && (si->early_stop)) break; // break from single domain case if rejected
         }
         i = -1;
         triggered = FALSE;
     }
   }
 
   /* Restore model to uni/multihit mode, and to its original length model */
   if (p7_IsMulti(save_mode)) p7_oprofile_ReconfigMultihit(om, saveL); 
   else                       p7_oprofile_ReconfigUnihit  (om, saveL); 
 
   return status;
 }

/**
 * @brief Determine if a region contains multiple domains based on occupancy thresholds.
 *
 * @param[in] ddef Domain definition object with occupancy stats.
 * @param[in]  i  Region start residue index.
 * @param[in]  j  Region end residue index.
 *
 * @return TRUE if expected domains >= rt3; FALSE otherwise.
 */
static int
is_multidomain_region(P7_DOMAINDEF *ddef, int i, int j)
{
  int   z;
  float max = -1.0;
  for (z = i; z <= j; z++) {
    float expected_n = ESL_MIN(ddef->etot[z] - ddef->etot[i-1],
                               ddef->btot[j] - ddef->btot[z-1]);
    max = ESL_MAX(max, expected_n);
  }
  return (max >= ddef->rt3) ? TRUE : FALSE;
}

/**
 * @brief Collect an ensemble of stochastic traces and cluster to define envelopes.
 *
 * @param[in]  ddef  Domain definition object.
 * @param[in]  om    Profile HMM.
 * @param[in]  dsq   Digitized sequence.
 * @param[in]  ireg  Region start.
 * @param[in]  jreg  Region end.
 * @param[in]  fwd   Forward matrix.
 * @param[in]  wrk   Workspace matrix.
 * @param[out] ret_nc Number of clusters returned.
 *
 * @return eslOK on success.
 */
static int
region_trace_ensemble(P7_DOMAINDEF *ddef, const P7_OPROFILE *om,
                      const ESL_DSQ *dsq, int ireg, int jreg,
                      const P7_OMX *fwd, P7_OMX *wrk, int *ret_nc)
{
   int    Lr  = jreg-ireg+1;
   int    t, d, d2;
   int    nov, n;
   int    nc;
   int    pos;
   float  null2[p7_MAXCODE];
 
   esl_vec_FSet(ddef->n2sc+ireg, Lr, 0.0); /* zero the null2 scores in region */
 
   /* By default, we make results reproducible by forcing a reset of
    * the RNG to its originally seeded state.
    */
   if (ddef->do_reseeding) 
     esl_randomness_Init(ddef->r, esl_randomness_GetSeed(ddef->r));
 
   /* Collect an ensemble of sampled traces; calculate null2 odds ratios from these */
   for (t = 0; t < ddef->nsamples; t++) {
       p7_StochasticTrace(ddef->r, dsq+ireg-1, Lr, om, fwd, ddef->tr);
       p7_trace_Index(ddef->tr);
 
       pos = 1;
       for (d = 0; d < ddef->tr->ndom; d++) {
         p7_spensemble_Add(ddef->sp, t, ddef->tr->sqfrom[d]+ireg-1, ddef->tr->sqto[d]+ireg-1, ddef->tr->hmmfrom[d], ddef->tr->hmmto[d]);
    
         p7_Null2_ByTrace(om, ddef->tr, ddef->tr->tfrom[d], ddef->tr->tto[d], wrk, null2);
        
         /* residues outside domains get bumped +1: because f'(x) = f(x), so f'(x)/f(x) = 1 in these segments */
         for (; pos <= ddef->tr->sqfrom[d]; pos++) ddef->n2sc[ireg+pos-1] += 1.0;
    
         /* Residues inside domains get bumped by their null2 ratio */
         for (; pos <= ddef->tr->sqto[d];   pos++) ddef->n2sc[ireg+pos-1] += null2[dsq[ireg+pos-1]];
      }

      /* the remaining residues in the region outside any domains get +1 */
      for (; pos <= Lr; pos++)  ddef->n2sc[ireg+pos-1] += 1.0;

      p7_trace_Reuse(ddef->tr);        
    }
 
   /* Convert the accumulated n2sc[] ratios in this region to log odds null2 scores on each residue. */
   for (pos = ireg; pos <= jreg; pos++)
     ddef->n2sc[pos] = logf(ddef->n2sc[pos] / (float) ddef->nsamples);
 
   /* Cluster the ensemble of traces to break region into envelopes. */
   p7_spensemble_Cluster(ddef->sp, ddef->min_overlap, ddef->of_smaller, ddef->max_diagdiff, ddef->min_posterior, ddef->min_endpointp, &nc);
 
   /* A little hacky now. Remove "dominated" domains relative to seq coords. */
   for (d = 0; d < nc; d++) 
     ddef->sp->assignment[d] = 0; /* overload <assignment> to flag that a domain is dominated */
 
   /* who dominates who? (by post prob) */
   for (d = 0; d < nc; d++) {
       for (d2 = d+1; d2 < nc; d2++) {
        nov = ESL_MIN(ddef->sp->sigc[d].j, ddef->sp->sigc[d2].j) - ESL_MAX(ddef->sp->sigc[d].i, ddef->sp->sigc[d2].i) + 1;
        if (nov == 0) break;
        n   = ESL_MIN(ddef->sp->sigc[d].j - ddef->sp->sigc[d].i + 1,  ddef->sp->sigc[d2].j - ddef->sp->sigc[d2].i + 1);
        if ((float) nov / (float) n >= 0.8) { /* overlap */
          if (ddef->sp->sigc[d].prob > ddef->sp->sigc[d2].prob) ddef->sp->assignment[d2] = 1;
          else                                                  ddef->sp->assignment[d]  = 1;
        }
      }
   }
       
   /* shrink the sigc list, removing dominated domains */
   d = 0;
   for (d2 = 0; d2 < nc; d2++) {
      if (ddef->sp->assignment[d2]) continue; /* skip domain d2, it's dominated. */
      if (d != d2) memcpy(ddef->sp->sigc + d, ddef->sp->sigc + d2, sizeof(struct p7_spcoord_s));
      d++;
    }
   ddef->sp->nc = d;
   *ret_nc = d;
   return eslOK;
 }

/**
 * @brief Rescore a single domain, calculate PID/E-value, and append rejected targets.
 *
 * @param[in]  ddef   Domain definition object.
 * @param[in]  om     Profile HMM.
 * @param[in]  sq     Target sequence.
 * @param[in]  qsq    Query sequence.
 * @param[in]  ox1    Forward scores matrix.
 * @param[in]  ox2    Backward/posterior matrix.
 * @param[in]  i      Domain start in target.
 * @param[in]  j      Domain end in target.
 * @param[in]  null2_is_done Flag if null2 already computed.
 * @param[in]  bg     Background model.
 * @param[in]  long_target Flag for long target.
 * @param[in]  bg_tmp Temp background.
 * @param[in]  scores_arr Scratch scoring array.
 * @param[in]  fwd_emissions_arr Scratch emissions.
 * @param[in]  result Result info struct.
 * @param[in]  nullsc Null model score ptr.
 * @param[in]  si    Global SLEDGE_INFO.
 *
 * @return eslOK on success, or error status.
 */
static int
rescore_domain_cust(P7_DOMAINDEF *ddef, P7_OPROFILE *om,
                    const ESL_SQ *sq, const ESL_SQ *qsq,
                    P7_OMX *ox1, P7_OMX *ox2,
                    int i, int j, int null2_is_done,
                    P7_BG *bg, int long_target,
                    P7_BG *bg_tmp,
                    float *scores_arr, float *fwd_emissions_arr,
                    RESULT_INFO *result, float *nullsc,
                    SLEDGE_INFO *si)
{
   P7_DOMAIN     *dom           = NULL;
   int            Ld            = j-i+1;
   float          domcorrection = 0.0;
   float          envsc, oasc;
   int            z;
   int            pos;
   float          null2[p7_MAXCODE];
   int            status;
   int            max_env_extra = 20;
   int            orig_L;
   float          eval;
   float          pid;
   char           decision;
 
   
   p7_Forward (sq->dsq + i-1, Ld, om, ox1, &envsc);
   p7_Backward(sq->dsq + i-1, Ld, om, ox1, ox2, NULL);
 
   status = p7_Decoding(om, ox1, ox2, ox2);      /* <ox2> is now overwritten with post probabilities     */
   if (status == eslERANGE) { /* rare: numeric overflow; domain is assumed to be repetitive garbage [J3/119-121] */
     goto ERROR;
   }
 
   /* Find an optimal accuracy alignment */
   p7_OptimalAccuracy(om, ox2, ox1, &oasc);      /* <ox1> is now overwritten with OA scores              */
   p7_OATrace        (om, ox2, ox1, ddef->tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */
 
   /* hack the trace's sq coords to be correct w.r.t. original dsq */
   for (z = 0; z < ddef->tr->N; z++)
     if (ddef->tr->i[z] > 0) ddef->tr->i[z] += i-1;
 
   /* Find the percent identity and accept/reject the sequence based on threshold */
   find_pid(ddef->tr, sq, qsq, &pid, &decision, si);
   
   /////////// IS THIS BOXED SECTION NEEDED? REMOVE, TEST AND UPDATE ///////////
   /* get ptr to next empty domain structure in domaindef's results */
   if (ddef->ndom == ddef->nalloc) {
     ESL_REALLOC(ddef->dcl, sizeof(P7_DOMAIN) * (ddef->nalloc*2));
     ddef->nalloc *= 2;
   }
   /////////////////////////////////////////////////////////////////////////////
 
  dom = &(ddef->dcl[ddef->ndom]);
  dom->ad             = p7_alidisplay_Create(ddef->tr, 0, om, sq, qsq);
  dom->scores_per_pos = NULL;

  if (!null2_is_done) {
    p7_Null2_ByExpectation(om, ox2, null2);
    for (pos = i; pos <= j; pos++)
      ddef->n2sc[pos]  = logf(null2[sq->dsq[pos]]);
  }
  for (pos = i; pos <= j; pos++)
    domcorrection += ddef->n2sc[pos]; /* domcorrection is in units of NATS */

  dom->domcorrection = domcorrection; /* in units of NATS */
  dom->ienv          = i;
  dom->jenv          = j;
  dom->envsc         = envsc;         /* in units of NATS */
  dom->oasc          = oasc;          /* in units of expected # of correctly aligned residues */
  dom->dombias       = 0.0;           /* gets set later, using bg->omega and dombias */
  dom->bitscore      = 0.0;           /* gets set later by caller, using envsc, null score, and dombias */
  dom->lnP           = 0.0;           /* gets set later by caller, using bitscore */
  dom->is_reported   = FALSE;         /* gets set later by caller */
  dom->is_included   = FALSE;         /* gets set later by caller */

  ddef->ndom++;
  
  /* E-value calculation */
  Ld = dom->jenv - dom->ienv + 1;
  dom->bitscore = dom->envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3));                 /* NATS, for the moment...   */
  dom->dombias  = p7_FLogsum(0.0, log(bg->omega) + dom->domcorrection);                             /* NATS, and will stay so    */
  dom->bitscore = (dom->bitscore - (*nullsc + dom->dombias)) / eslCONST_LOG2;                       /* now BITS, as it should be */
  dom->lnP      = esl_exp_logsurv (dom->bitscore,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);  /* p-value in log space      */
  eval = dom->lnP + log(si->db_size);                                                               /* E-value based on db size  */

  if ((eval < log(si->E)) && (decision == REJECT))                                                  /* E-value threshold check   */
    append_target(result, sq->name, pid, eval, ddef->tr->N);
  
   p7_trace_Reuse(ddef->tr);
 
   return eslOK;
 
  ERROR:
   p7_trace_Reuse(ddef->tr);
   return status;
 }

/**
 * @brief Compute percent identity from a traceback and decide accept/reject.
 *
 * @param[in]  tr       Trace alignment object.
 * @param[in]  sq       Target sequence.
 * @param[in]  qsq      Query sequence.
 * @param[out] pid      Computed percent identity.
 * @param[out] decision Character flag ACCEPT or REJECT.
 * @param[in]  si       Global SLEDGE_INFO with pid thresholds.
 *
 * @return eslOK on success.
 */
static int
find_pid(P7_TRACE *tr, const ESL_SQ *sq,
         const ESL_SQ *qsq, float *pid, char *decision,
         SLEDGE_INFO *si)
{
  int    j;
  double match = 0.0, q_length = 0.0, t_length = 0.0;

  /* Go through the trace */
  for (j = 0; j < tr->N; j++) {
    /* Find exact matches */
    if (tr->st[j] < p7T_M || tr->st[j] > p7T_I) continue;
    if (tr->st[j] == p7T_M && qsq->dsq[tr->k[j]] == sq->dsq[tr->i[j]]
        && qsq->dsq[tr->k[j]] < 20 && sq->dsq[tr->i[j]] < 20) {
      match++; q_length++;
    } else if (tr->st[j] == p7T_D) {
      t_length++;
    } else {
      q_length++;
    }
  }

  /* Compute pid and decision */
  *pid = match / q_length;
  *decision = (*pid > si->pid_low && *pid <= si->pid_high)
                ? ACCEPT : REJECT;

  return eslOK;
}