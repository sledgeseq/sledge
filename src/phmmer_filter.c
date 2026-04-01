/* phmmer: search a protein sequence against a protein database
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include <unistd.h>
#include "esl_threads.h"

#include "hmmer.h"
#include "sledge_dev.h"

/**
 * @name Command Line Options
 *
 * Defines all supported options for controlling output, scoring, thresholds,
 * acceleration heuristics, Sledgehmmer-specific behavior, and parallelism.
 */
static ESL_OPTIONS options[] = {
  /* name             type              default   env  range   toggles   reqs   incomp                             help                                       docgroup*/
  { "-h",             eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  NULL,              "show brief help on version and usage",                         1 },
/* Control of output */
  { "-A",             eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save multiple alignment of hits to file <f>",                  2 },
  { "--tblout",       eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",    eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout",   eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",          eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  NULL,              "prefer accessions over names in output",                       2 },
  { "--noali",        eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  NULL,              "don't output alignments, so output is smaller",                2 },
  { "--notextw",      eslARG_NONE,         NULL, NULL, NULL,        NULL,  NULL, "--textw",          "unlimit ASCII text output line width",                         2 },
  { "--textw",        eslARG_INT,         "120", NULL, "n>=120",    NULL,  NULL, "--notextw",        "set max width of ASCII text output lines",                     2 },
/* Control of scoring system */
  { "--popen",        eslARG_REAL,       "0.02", NULL, "0<=x<0.5",  NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",      eslARG_REAL,        "0.4", NULL, "0<=x<1",    NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mx",           eslARG_STRING, "BLOSUM62", NULL, NULL,        NULL,  NULL,  "--mxfile",        "substitution score matrix choice (of some built-in matrices)", 3 },
  { "--mxfile",       eslARG_INFILE,       NULL, NULL, NULL,        NULL,  NULL,  "--mx",            "read substitution score matrix from file <f>",                 3 },
/* Control of reporting thresholds */
  { "-E",             eslARG_REAL,       "10.0", NULL, "x>0",       NULL,  NULL,  REPOPTS,           "report sequences <= this E-value threshold in output",         4 },
  { "-T",             eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  REPOPTS,           "report sequences >= this score threshold in output",           4 },
  { "--domE",         eslARG_REAL,       "10.0", NULL, "x>0",       NULL,  NULL,  DOMREPOPTS,        "report domains <= this E-value threshold in output",           4 },
  { "--domT",         eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  DOMREPOPTS,        "report domains >= this score cutoff in output",                4 },
/* Control of inclusion thresholds */
  { "--incE",         eslARG_REAL,       "0.01", NULL, "x>0",       NULL,  NULL,  INCOPTS,           "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",         eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  INCOPTS,           "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",      eslARG_REAL,       "0.01", NULL, "x>0",       NULL,  NULL,  INCDOMOPTS,        "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",      eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  INCDOMOPTS,        "consider domains >= this score threshold as significant",      5 },
/* Model-specific thresholding for both reporting and inclusion (unused in phmmer)*/
  { "--cut_ga",       eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  THRESHOPTS,        "use profile's GA gathering cutoffs to set all thresholding",  99 },
  { "--cut_nc",       eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  THRESHOPTS,        "use profile's NC noise cutoffs to set all thresholding",      99 },
  { "--cut_tc",       eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  THRESHOPTS,        "use profile's TC trusted cutoffs to set all thresholding",    99 },
/* Control of filter pipeline */
  { "--max",          eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL, "--F1,--F2,--F3",   "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",           eslARG_REAL,       "0.02", NULL, NULL,        NULL,  NULL, "--max",            "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",           eslARG_REAL,       "1e-3", NULL, NULL,        NULL,  NULL, "--max",            "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",           eslARG_REAL,       "1e-5", NULL, NULL,        NULL,  NULL, "--max",            "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",       eslARG_NONE,        NULL,  NULL, NULL,        NULL,  NULL, "--max",            "turn off composition bias filter",                             7 },
/* Control of E-value calibration */
  { "--EmL",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",          eslARG_INT,         "100", NULL,"n>0",        NULL,  NULL,  NULL,              "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",          eslARG_REAL,       "0.04", NULL,"0<x<1",      NULL,  NULL,  NULL,              "tail mass for Forward exponential tail tau fit",              11 },   
/* Sledgehmmer specific options */  
  { "--qblock",       eslARG_INT,        "10",    NULL, NULL,       NULL,  NULL,  NULL,              "max queries handled at a time",                               12 },
  { "--tblock",       eslARG_INT,        "1000",  NULL, NULL,       NULL,  NULL,  NULL,              "max targets handled at a time",                               12 },
  { "--all_hits",   eslARG_NONE,         FALSE,   NULL, NULL,       NULL,  NULL,  NULL,              "early stopping toggle",                                       12 },
  { "--format",       eslARG_INT,         "1",    NULL, NULL,       NULL,  NULL,  NULL,              "output format",                                               12 },
  { "--halt",         eslARG_INT,        "-1",    NULL, NULL,       NULL,  NULL,  NULL,              "stop after n sequences - for debugging",                      12 },
  { "--suppress",     eslARG_NONE,       FALSE,   NULL, NULL,       NULL,  NULL,  NULL,              "turn off progress bar",                                       12 },
  { "--task_id",      eslARG_INT,        "0",     NULL, NULL,       NULL,  NULL,  NULL,              "slurm array task ID or shard number",                         12 },
  { "--shard_id",     eslARG_STRING,     "",      NULL, NULL,       NULL,  NULL,  NULL,              "second level shard number as a string",                       12 },
  { "--plow",         eslARG_REAL,       "0.00",  NULL, NULL,       NULL,  NULL,  NULL,              "pid lower limit for accepting a sequence",                    12 },
  { "--phigh",        eslARG_REAL,       "0.00",  NULL, NULL,       NULL,  NULL,  NULL,              "pid upper limit for accepting a sequence",                    12 },
  { "--qsize",        eslARG_INT,        "1",     NULL, NULL,       NULL,  NULL,  NULL,              "num queries per thread",                                      12 },
  { "-o",             eslARG_STRING,     "-",     NULL, NULL,       NULL,  NULL,  NULL,              "output file name",                                            12 },
  { "--freq",         eslARG_INT,        "1000",  NULL, NULL,       NULL,  NULL,  NULL,              "frequency of progress updates",                               12 },
  { "--cpu",          eslARG_INT,        "1","HMMER_NCPU", "n>=0",  NULL,  NULL,  CPUOPTS,           "number of CPU cores to use",                                  12 },
/* other options */
  { "--nonull2",      eslARG_NONE,        NULL,  NULL, NULL,        NULL,  NULL,  NULL,              "turn off biased composition score corrections",               12 },
  { "-Z",             eslARG_REAL,       FALSE,  NULL, "x>0",       NULL,  NULL,  NULL,              "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",         eslARG_REAL,       FALSE,  NULL, "x>0",       NULL,  NULL,  NULL,              "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",         eslARG_INT,         "42",  NULL, "n>=0",      NULL,  NULL,  NULL,              "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--qformat",      eslARG_STRING,      NULL,  NULL, NULL,        NULL,  NULL,  NULL,              "assert query <seqfile> is in format <s>: no autodetection",   12 },
  { "--tformat",      eslARG_STRING,      NULL,  NULL, NULL,        NULL,  NULL,  NULL,              "assert target <seqdb> is in format <s>>: no autodetection",   12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <qdb> <tdb>";
static char banner[] = "filter hits from qdb to tdb";

/**
 * @brief Parse and verify command-line arguments, display help if requested.
 *
 * This function initializes an ESL_GETOPTS object with the static
 * options[] array, processes the environment and argv, verifies the
 * configuration, and handles the '-h' help flag by printing usage
 * and grouped option descriptions.
 *
 * @param[in]  argc     Argument count (from main).
 * @param[in]  argv     Argument vector (from main).
 * @param[out] ret_go   Pointer to store the allocated ESL_GETOPTS object.
 * @param[out] ret_qfile Placeholder for query file argument (unused here).
 * @param[out] ret_dbfile Pointer to store the target database filename.
 *
 * @return eslOK on success; exits on help or fatal parse errors.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_qfile, char **ret_dbfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
  
      sledge_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      if (puts("\nSledgeHMMER and expert options:")                          < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80);
      
      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling scoring system:")                      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOptions controlling E value calibration:")                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);  
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)    { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_qfile  = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <qdb> argument on command line")     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_dbfile = esl_opt_GetArg(go, 2)) == NULL) { if (puts("Failed to get <tdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_qfile, "-") == 0 && strcmp(*ret_dbfile, "-") == 0) 
    { if (puts("Either <qdb> or <tdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");  goto FAILURE; }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere basic options are:")                                       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);

  exit(1);

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

/**
 * @brief Main program: orchestrate database filtering workflow.
 *
 * Reads the query and target databases then iteratively checks
 * for sequence pairs that fail the specified similarity criteria.
 * Progress is reported periodically, and output is given in a full
 * tabular form (considering all failed comparisons) or a list of
 * accept/reject IDs (considering only the first failed comparison).
 *
 * @param[in] argc Argument count.
 * @param[in] argv Argument vector.
 *
 * @return exit status code (0 on success).
 */
int
main(int argc, char **argv)
{ 
  int status = eslOK;

  ESL_GETOPTS     *go               = NULL;	              /* command line processing    */
  ESL_SQFILE      *dbfp             = NULL;               /* open dbfile                */
  ESL_SQFILE      *qfp              = NULL;               /* open qfile                 */
  ESL_SQ_BLOCK    *qdb              = NULL;               /* query database             */
  ESL_SQ_BLOCK    *tdb              = NULL;               /* target database            */
  ESL_SQ_BLOCK    *cdb              = NULL;               /* chunk database             */
  int              qformat          = eslSQFILE_UNKNOWN;  /* format of qfile            */
  int              tformat          = eslSQFILE_UNKNOWN;  /* format of tfile            */
  ESL_ALPHABET    *abc              = NULL;               /* sequence alphabet          */
  char            *results          = NULL;               /* results of seq inclusion   */
  struct cfg_s     cfg;                                   /* configuration data         */
  SLEDGE_INFO      si;                                    /* Sledgehmmer info           */
  int              num_records;                           /* Total sequences in qdb     */
  int              seq_ctr;                               /* Count sequences processed  */
  char            *out_file;                              /* Output file name           */
  int              len;                                   /* Calculate file name length */
  clock_t          start_time;                            /* Keep track of start time   */
  clock_t          elapsed_time;                          /* Elapsed time               */
  bool             suppress;                              /* Suppress progress update   */
  int              freq;                                  /* Freq of progress update    */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.qfile        = NULL;
  cfg.dbfile       = NULL;
  cfg.do_mpi       = FALSE;	         /* this gets reset below, if we init MPI */
  cfg.nproc        = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank      = 0;		           /* this gets reset below, if we init MPI */
  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;

  /* Initializations */

  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &cfg.qfile, &cfg.dbfile);

  /* Initialize parameters for Sledgehmmer*/
  si.pid_low      = esl_opt_GetReal    (go, "--plow");
  si.pid_high     = esl_opt_GetReal    (go, "--phigh");
  si.qsize        = esl_opt_GetInteger (go, "--qsize");
  si.cores        = esl_opt_GetInteger (go, "--cpu");
  si.E            = esl_opt_GetReal    (go, "-E");
  si.db_size      = esl_opt_GetReal    (go, "-Z");
  si.out_path     = esl_opt_GetString  (go, "-o");
  si.task_id      = esl_opt_GetInteger (go, "--task_id");
  si.format       = esl_opt_GetInteger (go, "--format");
  suppress        = esl_opt_GetBoolean (go, "--suppress");
  si.offset       = 0;
  si.out_fp       = NULL;
  
  if(esl_opt_GetBoolean (go, "--all_hits"))
    si.early_stop = FALSE;
  else
    si.early_stop = TRUE;

  /* Set unused Sledgehmmer parameters to NULL */
  si.train_only   = FALSE;
  si.test_only    = FALSE;
  si.train_frac   = 0.0;
  si.init_chunk   = 0;
  si.test_limit   = 0;
  si.train_path   = NULL;
  si.test_path    = NULL;
  si.discard_path = NULL;

  /* If caller declared input formats, decode them */
  if (esl_opt_IsOn(go, "--qformat")){
    qformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (qformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }
  if (esl_opt_IsOn(go, "--tformat")) {
    tformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (tformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }
  
  /* Initialize alphabet */
  abc = esl_alphabet_Create(eslAMINO);

  /* Open and read the target database file */
  status =  esl_sqfile_OpenDigital(abc, cfg.dbfile, tformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open target sequence database %s for reading\n",      cfg.dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Target sequence database file %s is empty or misformatted\n",   cfg.dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg.dbfile);
  
  tdb = esl_sq_CreateDigitalBlock(esl_opt_GetInteger(go, "--tblock"), abc);
  if (tdb == NULL) p7_Fail("Failed to allocate sequence database block");

  status = read_cust(dbfp, tdb);
  if (status == eslENOSPACE)
    p7_Fail("Failed to read target sequences: --tblock %d is smaller than the number of sequences in the database\n",
            esl_opt_GetInteger(go, "--tblock"));
  if (status != eslOK) p7_Fail("Failed to read target sequences from database");
  esl_sqfile_Close(dbfp);

  /* Open and read the query database file */
  status =  esl_sqfile_OpenDigital(abc, cfg.qfile, qformat, p7_SEQDBENV, &qfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open query sequence database %s for reading\n",       cfg.qfile);
  else if (status == eslEFORMAT)   p7_Fail("Query sequence database file %s is empty or misformatted\n",    cfg.qfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg.qfile);

  qdb = esl_sq_CreateDigitalBlock(esl_opt_GetInteger(go, "--qblock"), abc);
  if (qdb == NULL) p7_Fail("Failed to allocate sequence database block");

  status = read_cust(qfp, qdb);
  if (status == eslENOSPACE)
    p7_Fail("Failed to read query sequences: --qblock %d is smaller than the number of sequences in the query file\n",
            esl_opt_GetInteger(go, "--qblock"));
  if (status != eslOK) p7_Fail("Failed to read query sequences from database");
  esl_sqfile_Close(qfp);

  /* Allocate a database to process each chunk of queries */
  cdb = esl_sq_CreateDigitalBlock(si.cores * si.qsize, abc); /* Allocate for max possible */
  if (cdb == NULL) p7_Fail("Failed to allocate chunk database block");

  /* Calculate total number of queries to process */
  num_records = qdb->count;

  /* Allocate large buffer for filter results */
  if (si.format == R_STRING) ESL_ALLOC(results, (num_records+1) * sizeof(char));
  else  ESL_ALLOC(results, BUFFER_MAX * sizeof(char));
  
  /* Open the output file (or stdout) for permanent writing */
  if (strcmp(si.out_path, "-") == 0)
    si.out_fp = stdout;
  else {
    len = snprintf(NULL, 0, "%s_%d%s.txt", si.out_path, si.task_id, esl_opt_GetString(go, "--shard_id")) + 1;
    ESL_ALLOC(out_file, len*sizeof(char));
    snprintf(out_file, len, "%s_%d%s.txt", si.out_path, si.task_id, esl_opt_GetString(go, "--shard_id"));
    if ((si.out_fp = fopen(out_file, "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n", out_file);
  }

  /* Header in output file/stream if printing tabular output */
  if (si.format == R_ALL) {
    fprintf(si.out_fp, "%-15s %-15s %6s %7s %12s %8s\n",
      "Query", "Target", "PID", "Length", "E-value", "Decision");

    fprintf(si.out_fp, "--------------- --------------- ------ ------- ------------ --------\n"); /* Just pretty printing */
  }

  start_time = get_coarse_time();
  freq = esl_opt_GetInteger (go, "--freq");

  /* Loop through all sequences in database */

  for (seq_ctr = 0; seq_ctr < num_records; seq_ctr++) {
    add_seq(cdb, qdb->list + seq_ctr, abc);

    /* Check of query size is reached or end of db is reached */
    if (((cdb->count == si.cores * si.qsize) || (seq_ctr == num_records-1))) {
      if (seq_ctr == num_records-1)
        si.cores = (cdb->count + si.qsize - 1) / si.qsize;

      status = assign_master(go, cdb, tdb, &si, results);
      if (si.format == R_STRING) fwrite(results, 1, cdb->count, si.out_fp);

      cdb->count = 0;
    }

    /* Display progress update based on frequency */
    if ((!suppress) && (seq_ctr % freq == 0)) {
      elapsed_time = get_coarse_time() - start_time;
      print_progress(seq_ctr+1, num_records, elapsed_time);

      fflush(stderr);
    }
  }

  /* Flush the buffer one last time if buffer isn't empty */
  if (si.offset > 0)
    fwrite(results, 1, si.offset, si.out_fp);

  /* Print progress for one final time */
  elapsed_time = get_coarse_time() - start_time;
  if (!suppress) print_progress(seq_ctr, num_records, elapsed_time);

  /* Close all files if used */

  if (si.out_fp != stdout) {
    fclose(si.out_fp);
    free(out_file);
  }

  /* Mop the floor and leave */
  
  free(results);
  
  sledge_worker_pool_destroy();

  esl_sq_DestroyBlock(qdb);
  esl_sq_DestroyBlock(tdb);
  esl_sq_DestroyBlock(cdb);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);

  ERROR:
    exit(status);

  exit(status);
}
