/**
 * @file sledge_splitter.c
 * @brief Split a protein sequence database into training, testing, and optional validation sets using a profmark-style algorithm.
 *
 * This tool reads a sequence database, partitions it according to random sampling
 * and configurable proportions, then writes out FASTA files for train, test,
 * (and optionally validation) sets. It supports iterative, chunked processing
 * with multithreaded phmmer-based assignment of sequences between sets.
 *
 * Usage: sledge_splitter [-options] <seqdb>
 *
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
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
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
  { "-h",             eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  NULL,               "show brief help on version and usage",                         1 },
/* Control of output */
  { "-A",             eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,               "save multiple alignment of hits to file <f>",                  2 },
  { "--tblout",       eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,               "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",    eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,               "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout",   eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,               "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",          eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  NULL,               "prefer accessions over names in output",                       2 },
  { "--noali",        eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  NULL,               "don't output alignments, so output is smaller",                2 },
  { "--notextw",      eslARG_NONE,         NULL, NULL, NULL,      NULL,  NULL, "--textw",           "unlimit ASCII text output line width",                         2 },
  { "--textw",        eslARG_INT,         "120", NULL, "n>=120",  NULL,  NULL, "--notextw",         "set max width of ASCII text output lines",                     2 },
/* Control of scoring system */
  { "--popen",        eslARG_REAL,       "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,               "gap open probability",                                         3 },
  { "--pextend",      eslARG_REAL,        "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,               "gap extend probability",                                       3 },
  { "--mx",           eslARG_STRING, "BLOSUM62", NULL, NULL,      NULL,  NULL,  "--mxfile",         "substitution score matrix choice (of some built-in matrices)", 3 },
  { "--mxfile",       eslARG_INFILE,       NULL, NULL, NULL,      NULL,  NULL,  "--mx",             "read substitution score matrix from file <f>",                 3 },
/* Control of reporting thresholds */
  { "-E",             eslARG_REAL,       "10.0", NULL, "x>0",     NULL,  NULL,  REPOPTS,            "report sequences <= this E-value threshold in output",         4 },
  { "-T",             eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  REPOPTS,            "report sequences >= this score threshold in output",           4 },
  { "--domE",         eslARG_REAL,       "10.0", NULL, "x>0",     NULL,  NULL,  DOMREPOPTS,         "report domains <= this E-value threshold in output",           4 },
  { "--domT",         eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  DOMREPOPTS,         "report domains >= this score cutoff in output",                4 },
/* Control of inclusion thresholds */
  { "--incE",         eslARG_REAL,       "0.01", NULL, "x>0",     NULL,  NULL,  INCOPTS,            "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",         eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  INCOPTS,            "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",      eslARG_REAL,       "0.01", NULL, "x>0",     NULL,  NULL,  INCDOMOPTS,         "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",      eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  INCDOMOPTS,         "consider domains >= this score threshold as significant",      5 },
/* Model-specific thresholding for both reporting and inclusion (unused in phmmer)*/
  { "--cut_ga",       eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  THRESHOPTS,         "use profile's GA gathering cutoffs to set all thresholding",  99 },
  { "--cut_nc",       eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  THRESHOPTS,         "use profile's NC noise cutoffs to set all thresholding",      99 },
  { "--cut_tc",       eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  THRESHOPTS,         "use profile's TC trusted cutoffs to set all thresholding",    99 },
/* Control of filter pipeline */
  { "--max",          eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3",    "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",           eslARG_REAL,       "0.02", NULL, NULL,      NULL,  NULL, "--max",             "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",           eslARG_REAL,       "1e-3", NULL, NULL,      NULL,  NULL, "--max",             "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",           eslARG_REAL,       "1e-5", NULL, NULL,      NULL,  NULL, "--max",             "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",       eslARG_NONE,        NULL,  NULL, NULL,      NULL,  NULL, "--max",             "turn off composition bias filter",                             7 },
/* Control of E-value calibration */
  { "--EmL",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,               "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,               "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,               "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,               "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",          eslARG_INT,         "100", NULL,"n>0",      NULL,  NULL,  NULL,               "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,               "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",          eslARG_REAL,       "0.04", NULL,"0<x<1",    NULL,  NULL,  NULL,               "tail mass for Forward exponential tail tau fit",              11 },   
/* Sledgehmmer specific options */
  { "--qblock",       eslARG_INT,        "10",    NULL, NULL,      NULL,  NULL,  NULL,              "max queries handled at a time",                               12 },
  { "--dbblock",      eslARG_INT,        "1000",  NULL, NULL,      NULL,  NULL,  NULL,              "max targets handled at a time",                               12 },
  { "--halt",         eslARG_INT,        "-1",    NULL, NULL,      NULL,  NULL,  NULL,              "stop after n sequences - for debugging",                      12 },
  { "--resume",       eslARG_INT,        "0",     NULL, NULL,      NULL,  NULL,  NULL,              "resume from nth sequence",                                    12 },
  { "--test_limit",   eslARG_INT,        "75000", NULL, NULL,      NULL,  NULL,  NULL,              "required minimum number of test sequences",                   12 },
  { "--val_limit",    eslARG_INT,        "10000", NULL, NULL,      NULL,  NULL,  NULL,              "required minimum number of validation sequences",             12 },
  { "--init_chunk",   eslARG_INT,        "10",    NULL, NULL,      NULL,  NULL,  NULL,              "query chunk size till test limit",                            12 },
  { "--suppress",     eslARG_NONE,       FALSE,   NULL, NULL,      NULL,  NULL,  NULL,              "turn off progress bar",                                       12 },
  { "--task_id",      eslARG_INT,        "0",     NULL, NULL,      NULL,  NULL,  NULL,              "slurm array task ID or shard number",                         12 },
  { "--shard_id",     eslARG_STRING,     "",      NULL, NULL,      NULL,  NULL,  NULL,              "second level shard number as a string",                       12 },
  { "--plow",         eslARG_REAL,       "0.00",  NULL, NULL,      NULL,  NULL,  NULL,              "pid lower limit for accepting a sequence",                    12 },
  { "--phigh",        eslARG_REAL,       "0.00",  NULL, NULL,      NULL,  NULL,  NULL,              "pid upper limit for accepting a sequence",                    12 },
  { "--output_dir",   eslARG_STRING,      NULL,   NULL, NULL,      NULL,  NULL,  NULL,              "write train/test/val/discard to <dir> with default names",    12 },
  { "--train_path",   eslARG_STRING,     "train", NULL, NULL,      NULL,  NULL,  NULL,              "output basename prefix for train",                            12 },
  { "--test_path",    eslARG_STRING,     "test",        NULL, NULL, NULL,  NULL,  NULL,             "output basename prefix for test",                             12 },
  { "--val_path",     eslARG_STRING,     "val",         NULL, NULL, NULL,  NULL,  NULL,             "output basename prefix for val",                              12 },
  { "--discard_path", eslARG_STRING,     "discard",     NULL, NULL, NULL,  NULL,  NULL,             "output basename prefix for discard",                          12 },
  { "--load_tr",      eslARG_STRING,     "-",    NULL, NULL,       NULL,  NULL,  NULL,              "Resume from this train db",                                   12 },
  { "--load_te",      eslARG_STRING,     "-",    NULL, NULL,       NULL,  NULL,  NULL,              "Resume from this test db",                                    12 },
  { "--load_val",     eslARG_STRING,     "-",    NULL, NULL,       NULL,  NULL,  NULL,              "Resume from this val db",                                     12 },
  { "--freq",         eslARG_INT,        "1000",  NULL, NULL,      NULL,  NULL,  NULL,              "frequency of progress updates",                               12 },
  { "--cpu",          eslARG_INT,        "1","HMMER_NCPU", "n>=0", NULL,  NULL,  CPUOPTS,           "number of CPU cores to use",                                  12 },
  { "-o",             eslARG_STRING,     "-",     NULL, NULL,      NULL,  NULL,  NULL,              "output file name",                                            12 },
  /* other options */
  { "--nonull2",      eslARG_NONE,        NULL,  NULL, NULL,      NULL,  NULL,  NULL,               "turn off biased composition score corrections",               12 },
  { "-Z",             eslARG_REAL,       FALSE,  NULL, "x>0",     NULL,  NULL,  NULL,               "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",         eslARG_REAL,       FALSE,  NULL, "x>0",     NULL,  NULL,  NULL,               "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",         eslARG_INT,         "42",  NULL, "n>=0",    NULL,  NULL,  NULL,               "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--qformat",      eslARG_STRING,      NULL,  NULL, NULL,      NULL,  NULL,  NULL,               "assert query <seqfile> is in format <s>: no autodetection",   12 },
  { "--tformat",      eslARG_STRING,      NULL,  NULL, NULL,      NULL,  NULL,  NULL,               "assert target <seqdb> is in format <s>>: no autodetection",   12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqdb>";
static char banner[] = "split a database into train, test and optional validation set";

/* Resolve train/test/val/discard/stats path prefixes: optional --output_dir with default
 * basenames, or explicit --train_path/--test_path/... overrides. When --output_dir is set,
 * -o (stats) is placed under that directory unless -o is "-" (stdout). */
static void
resolve_output_paths(ESL_GETOPTS *go, SLEDGE_INFO *si)
{
  static char storage[5][4096];
  const char *dir = NULL;
  const char *o;

  if (esl_opt_IsOn(go, "--output_dir"))
    dir = esl_opt_GetString(go, "--output_dir");

  if (dir != NULL && dir[0] != '\0') {
    if (esl_opt_IsDefault(go, "--train_path"))
      snprintf(storage[0], sizeof(storage[0]), "%s/train", dir);
    else
      snprintf(storage[0], sizeof(storage[0]), "%s", esl_opt_GetString(go, "--train_path"));
    si->train_path = storage[0];

    if (esl_opt_IsDefault(go, "--test_path"))
      snprintf(storage[1], sizeof(storage[1]), "%s/test", dir);
    else
      snprintf(storage[1], sizeof(storage[1]), "%s", esl_opt_GetString(go, "--test_path"));
    si->test_path = storage[1];

    if (esl_opt_IsDefault(go, "--val_path"))
      snprintf(storage[2], sizeof(storage[2]), "%s/val", dir);
    else
      snprintf(storage[2], sizeof(storage[2]), "%s", esl_opt_GetString(go, "--val_path"));
    si->val_path = storage[2];

    if (esl_opt_IsDefault(go, "--discard_path"))
      snprintf(storage[3], sizeof(storage[3]), "%s/discard", dir);
    else
      snprintf(storage[3], sizeof(storage[3]), "%s", esl_opt_GetString(go, "--discard_path"));
    si->discard_path = storage[3];

    o = esl_opt_GetString(go, "-o");
    if (o != NULL && strcmp(o, "-") != 0) {
      /* Keep stats under output_dir: absolute -o → basename only; relative paths with ".."
       * → basename only so "../x" cannot escape the directory. */
      int strip_to_basename = (o[0] == '/' || strstr(o, "..") != NULL);
      if (strip_to_basename) {
        const char *slash = strrchr(o, '/');
        const char *name  = (slash && slash[1] != '\0') ? slash + 1 : o;
        if (name[0] == '\0' || strcmp(name, "..") == 0)
          name = "stats";
        snprintf(storage[4], sizeof(storage[4]), "%s/%s", dir, name);
      } else
        snprintf(storage[4], sizeof(storage[4]), "%s/%s", dir, o);
      si->out_path = storage[4];
    } else
      si->out_path = (char *)((o != NULL) ? o : "-");
  } else {
    si->train_path   = esl_opt_GetString(go, "--train_path");
    si->test_path    = esl_opt_GetString(go, "--test_path");
    si->val_path     = esl_opt_GetString(go, "--val_path");
    si->discard_path = esl_opt_GetString(go, "--discard_path");
    si->out_path     = esl_opt_GetString(go, "-o");
  }
}

/* Create --output_dir if missing. If the path exists, it must already be a directory (never replaced). */
static int
ensure_output_dir(const char *path)
{
  struct stat st;

  if (path == NULL || path[0] == '\0')
    return eslOK;
  if (stat(path, &st) == 0) {
    if (!S_ISDIR(st.st_mode))
      p7_Fail("--output_dir %s exists but is not a directory\n", path);
    return eslOK;
  }
  if (errno != ENOENT)
    p7_Fail("cannot stat --output_dir %s: %s\n", path, strerror(errno));
  if (mkdir(path, 0755) != 0)
    p7_Fail("Failed to create output directory %s: %s\n", path, strerror(errno));
  return eslOK;
}

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

      if (puts("\nSledgehmmer and expert options:")                          < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
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

  if (esl_opt_ArgNumber(go)                 != 1)    { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_dbfile = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_dbfile, "-") == 0) 
    { if (puts("Either <seqfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");  goto FAILURE; }

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
 * @brief Main program: orchestrate database splitting workflow.
 *
 * Reads the input sequence database, optionally loads existing train/test
 * blocks, initializes random sampling parameters, then iteratively
 * partitions the database in chunks. At each chunk, sequences are
 * assigned to train/test(/val) using phmmer-based assignment or random
 * sampling. Progress is reported periodically, and final FASTA files
 * are written for each partition.
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
  ESL_SQ_BLOCK    *train_db         = NULL;               /* train database             */
  ESL_SQ_BLOCK    *test_db          = NULL;               /* test database              */
  ESL_SQ_BLOCK    *val_db           = NULL;               /* val database               */
  ESL_SQ_BLOCK    *qdb              = NULL;               /* query database             */
  ESL_SQ_BLOCK    *cdb              = NULL;               /* two-group candidate db     */
  ESL_SQ_BLOCK    *discard_db       = NULL;               /* batched discard writes     */
  ESL_SQ_BLOCK    *sqdb             = NULL;               /* sequence database to split */
  int              dbformat         = eslSQFILE_UNKNOWN;  /* format of dbfile           */
  ESL_ALPHABET    *abc              = NULL;               /* sequence alphabet          */
  char            *results          = NULL;               /* results of seq inclusion   */
  struct cfg_s     cfg;                                   /* configuration data         */
  SLEDGE_INFO      si;                                    /* Sledgehmmer info           */
  int              num_records;                           /* Total sequences in db      */
  ESL_RANDOMNESS  *rnd;                                   /* Random number generator    */
  float            rnd_value;                             /* Random number value        */
  int              seq_ctr;                               /* Count sequences processed  */
  int              tmp_ctr;                               /* Temp counter               */
  int              resume;                                /* Resume from seq            */
  bool             assign_train;                          /* Flag for train assignment  */
  bool             write_train;                           /* Flag to save train seqs    */
  bool             write_test;                            /* Flag to save test seqs     */
  bool             write_val;                             /* Flag to save val seqs      */
  bool             test_done;                             /* Test set completed flag    */
  bool             val_done;                              /* Val set completed flag     */
  bool             suppress;                              /* Flag to suppress outputs   */
  int              Q_SIZE;                                /* Query chunk size           */
  double           current_assign_frac;                   /* Current train proportion   */
  char            *out_file;                              /* Output file name           */
  char            *discard_file;                          /* Discard sequence file name */
  FILE            *discard_fp;                            /* Discard sequence file ptr  */
  int              len;                                   /* Calculate file name length */
  double           start_time;                            /* Keep track of start time   */
  double           elapsed_time;                          /* Elapsed time               */
  int              freq;                                  /* Freq of progress update    */
  int              eff_test_limit;                        /* test_limit, offset if val merged into test */
  int              eff_val_limit;                         /* val_limit, offset if test merged into val  */

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

  /* Initialize parameters for Sledgehmmer */
  si.pid_low      = esl_opt_GetReal    (go, "--plow");
  si.pid_high     = esl_opt_GetReal    (go, "--phigh");
  si.cores        = esl_opt_GetInteger (go, "--cpu");
  si.E            = esl_opt_GetReal    (go, "-E");
  si.init_chunk   = esl_opt_GetInteger (go, "--init_chunk");
  si.test_limit   = esl_opt_GetInteger (go, "--test_limit");
  si.val_limit    = esl_opt_GetInteger (go, "--val_limit");
  si.task_id      = esl_opt_GetInteger (go, "--task_id");
  suppress        = esl_opt_GetBoolean (go, "--suppress");
  /* -Z default in options[] is FALSE (NULL): esl_opt_GetReal would read atof(NULL). When -Z is
   * unset (IsOn false), match E-value database size to --dbblock; otherwise use the given -Z. */
  if (! esl_opt_IsOn(go, "-Z"))
    si.db_size = esl_opt_GetInteger(go, "--dbblock");
  else
    si.db_size = (int) esl_opt_GetReal(go, "-Z");
  resolve_output_paths(go, &si);
  if (esl_opt_IsOn(go, "--output_dir"))
    ensure_output_dir(esl_opt_GetString(go, "--output_dir"));
  si.out_fp       = NULL;
  si.early_stop   = TRUE;  /* The splitter does early stopping             */
  si.format       = 0;     /* Store ACCEPT/REJECT string                   */
  si.train_skip_lo = -1;
  si.train_skip_hi = -1;
  si.merged_skip_lo = -1;
  si.merged_skip_hi = -1;
  si.merged_skip_in_test = 0;

  /* If caller declared input formats, decode them */
  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }
  
  /* Initialize alphabet */
  abc = esl_alphabet_Create(eslAMINO);
  
  /* Open the sequence database file to split */
  status =  esl_sqfile_OpenDigital(abc, cfg.dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open target sequence database %s for reading\n",      cfg.dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Target sequence database file %s is empty or misformatted\n",   cfg.dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg.dbfile);

  /* Create database block */
  sqdb = esl_sq_CreateDigitalBlock(esl_opt_GetInteger(go, "--dbblock"), abc);
  if (sqdb == NULL) p7_Fail("Failed to allocate sequence database block");

  /* Read sequences into database */
  status = read_cust(dbfp, sqdb);
  if (status != eslOK) p7_Fail("Failed to read target sequences from database");
  esl_sqfile_Close(dbfp);
  
  /* Calculate total number of sequences to process */
  num_records = sqdb->count;
  if (esl_opt_GetInteger(go, "--halt") > 0) num_records = esl_opt_GetInteger(go, "--halt");

  /* Allocate memory for the train SQ_BLOCK */
  train_db = esl_sq_CreateDigitalBlock((int) (5 * (si.test_limit+si.val_limit)), abc);
  if (train_db == NULL) p7_Fail("Failed to allocate train database block");
  

  /* Check if we are loading an existing train db initially */
  if (strcmp(esl_opt_GetString(go, "--load_tr"), "-") != 0) {
    char *fname = NULL;
    fname = esl_opt_GetString(go, "--load_tr");
    ESL_SQFILE *fp = NULL;
    status =  esl_sqfile_OpenDigital(abc, fname, dbformat, p7_SEQDBENV, &fp);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open train sequence database %s for reading\n",      fname);
    else if (status == eslEFORMAT)   p7_Fail("Train sequence database file %s is empty or misformatted\n",   fname);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening train sequence database file %s\n", status, fname);

    read_cust(fp, train_db);
    if (status != eslOK) p7_Fail("Failed to read train sequences from database");
    esl_sqfile_Close(fp);
  }

  /* Allocate memory for the test SQ_BLOCK */
  test_db = esl_sq_CreateDigitalBlock((int) (si.val_limit+si.init_chunk), abc);
  if (test_db == NULL) p7_Fail("Failed to allocate test database block");
  
  /* Check if we are loading an existing test db initially */
  if (strcmp(esl_opt_GetString(go, "--load_te"), "-") != 0) {
    char *fname = NULL;
    fname = esl_opt_GetString(go, "--load_te");
    ESL_SQFILE *fp = NULL;
    status =  esl_sqfile_OpenDigital(abc, fname, dbformat, p7_SEQDBENV, &fp);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open test sequence database %s for reading\n",      fname);
    else if (status == eslEFORMAT)   p7_Fail("Test sequence database file %s is empty or misformatted\n",   fname);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening test sequence database file %s\n", status, fname);
    
    read_cust(fp, test_db);
    if (status != eslOK) p7_Fail("Failed to read test sequences from database");
    esl_sqfile_Close(fp);
  }

  /* Allocate memory for the val SQ_BLOCK */
  val_db = esl_sq_CreateDigitalBlock((int) (si.val_limit+si.init_chunk), abc);
  if (val_db == NULL) p7_Fail("Failed to allocate val database block");

  /* Check if we are loading an existing val db initially */
  if (strcmp(esl_opt_GetString(go, "--load_val"), "-") != 0) {
    char *fname = NULL;
    fname = esl_opt_GetString(go, "--load_val");
    ESL_SQFILE *fp = NULL;
    status =  esl_sqfile_OpenDigital(abc, fname, dbformat, p7_SEQDBENV, &fp);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open val sequence database %s for reading\n",      fname);
    else if (status == eslEFORMAT)   p7_Fail("Val sequence database file %s is empty or misformatted\n",   fname);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening val sequence database file %s\n", status, fname);

    read_cust(fp, val_db);
    if (status != eslOK) p7_Fail("Failed to read val sequences from database");
    esl_sqfile_Close(fp);
  }

  /* Set initial cores and query chunk sizes */
  Q_SIZE = esl_opt_GetInteger(go, "--init_chunk");  /* Initial query set size        */
  current_assign_frac = 0.33;  /* 33-33-33 triple phase until a limit is reached */

  /* Allocate memory for the query SQ_BLOCK */
  qdb = esl_sq_CreateDigitalBlock(Q_SIZE, abc); /* Allocate for max possible, not just Q_SIZE */
  if (qdb == NULL) p7_Fail("Failed to allocate query database block");

  write_val = TRUE;
  write_train = TRUE;
  write_test = TRUE;
  eff_test_limit = si.test_limit;
  eff_val_limit  = si.val_limit;
  test_done = (test_db->count >= eff_test_limit);
  val_done  = (val_db->count >= eff_val_limit);
  if (val_done) current_assign_frac = 0.5;

  /* If a preloaded set already meets a limit, persist canonical output immediately. */
  if (test_done && write_test) {
    int skl = -1, skh = -1;
    if (si.merged_skip_lo >= 0 && si.merged_skip_in_test) {
      skl = si.merged_skip_lo;
      skh = si.merged_skip_hi;
    }
    save_seqs(NULL, si.test_path, si.task_id, test_db, skl, skh);
    write_test = FALSE;
  }
  if (val_done && write_val) {
    int skl = -1, skh = -1;
    if (si.merged_skip_lo >= 0 && !si.merged_skip_in_test) {
      skl = si.merged_skip_lo;
      skh = si.merged_skip_hi;
    }
    save_seqs(NULL, si.val_path, si.task_id, val_db, skl, skh);
    write_val = FALSE;
  }

  if (si.cores < Q_SIZE) {
    if (Q_SIZE % si.cores == 0) {
      si.qsize = (float) Q_SIZE / si.cores;
    }

    else {
      printf("Since cores < query chunk size, they must be exactly divisible\n");
      exit(-20);
    }
  }

  else {
    si.qsize = 1;       /* Feasible, so keep qsize 1 and cores as needed */
    si.cores = Q_SIZE;  /* Cores as much as needed for Q_SIZE                  */
  }

  /* Allocate max possible memory for results + 1 for \0 */
  ESL_ALLOC(results, (qdb->listSize+1)*sizeof(char));
  cdb = esl_sq_CreateDigitalBlock(qdb->listSize, abc);
  discard_db = esl_sq_CreateDigitalBlock(qdb->listSize, abc);
  if (cdb == NULL || discard_db == NULL)
    p7_Fail("Failed to allocate temporary sequence buffers");

  /* Create seeded random generator for the rest of the code to use */
  rnd = esl_randomness_Create(esl_opt_GetInteger (go, "--seed"));

  /* Open the discard file for permanent writing */

  len = snprintf(NULL, 0, "%s_%d.fasta", si.discard_path, si.task_id) + 1;
  ESL_ALLOC(discard_file, len*sizeof(char));
  
  snprintf(discard_file, len, "%s_%d.fasta", si.discard_path, si.task_id);
  if ((discard_fp = fopen(discard_file, "w")) == NULL) p7_Fail("Failed to open discard file %s for writing\n", discard_file);

  start_time = get_coarse_time();
  freq       = esl_opt_GetInteger (go, "--freq");
  resume     = esl_opt_GetInteger (go, "--resume");

  /* Loop through all sequences in database */
  for (seq_ctr = 0; seq_ctr < num_records; seq_ctr++) {
    if (test_done && val_done) break;

    if (seq_ctr < resume) continue; /* Resume code from said sequence number in db */
    
    else {
      add_seq(qdb, sqdb->list + seq_ctr, abc);
      
      /* Process queries in normal mode */
      if ((qdb->count== Q_SIZE) || (seq_ctr == num_records-1)) {
        rnd_value = esl_random(rnd);

        if (!test_done && !val_done) {
          if (rnd_value <= current_assign_frac)
            check_triple(qdb, train_db, test_db, val_db, discard_fp, &si, abc, go);
          else if (rnd_value > 2*current_assign_frac)
            check_triple(qdb, val_db, train_db, test_db, discard_fp, &si, abc, go);
          else
            check_triple(qdb, test_db, val_db, train_db, discard_fp, &si, abc, go);
        } else if (!test_done && val_done) {
          /* test_db and train_db already include appended val after val limit (see housekeeping). */
          assign_train = (rnd_value < 0.5);
          status = process_two_group_chunk(go, qdb, train_db, test_db, test_db, train_db,
                                           cdb, discard_db, discard_fp, abc, &si, assign_train, results);
        } else if (test_done && !val_done) {
          /* val_db and train_db already include appended test if test finished first. */
          assign_train = (rnd_value < 0.5);
          status = process_two_group_chunk(go, qdb, train_db, val_db, val_db, train_db,
                                           cdb, discard_db, discard_fp, abc, &si, assign_train, results);
        }
      }
      
      /* Housekeeping after processing query set */
      if (qdb->count == Q_SIZE) {
        qdb->count = 0;

        if (num_records - seq_ctr - 1 < Q_SIZE) {
          Q_SIZE = num_records - seq_ctr - 1;
          if (Q_SIZE <= 0) continue;

          if (si.cores >= Q_SIZE) {
            si.qsize = 1;                         /* Last chunk, use only required qsize           */
            si.cores = Q_SIZE;                    /* Last chunk, use only required cores           */
          }

          else {
            si.qsize = ceil((float) Q_SIZE / si.cores);   /* Core limited, so choose qsize wisely          */
            si.cores = ceil((float) Q_SIZE / si.qsize);   /* New cores based on adjusted qsize             */
          }
        }
      }

      /* Housekeeping after reaching val limit */
      if ((!val_done) && (val_db->count >= eff_val_limit)) {
        if (write_val) {
          int skl = -1, skh = -1;
          if (si.merged_skip_lo >= 0 && !si.merged_skip_in_test) {
            skl = si.merged_skip_lo;
            skh = si.merged_skip_hi;
          }
          save_seqs(NULL, si.val_path, si.task_id, val_db, skl, skh);
          write_val = FALSE;
        }
        /* Val reached before test: embed val into train and test, then search uses those blocks. */
        if (!test_done) {
          si.train_skip_lo = train_db->count;
          for (tmp_ctr = 0; tmp_ctr < val_db->count; tmp_ctr++)
            add_seq(train_db, val_db->list + tmp_ctr, abc);
          si.train_skip_hi = train_db->count;
          si.merged_skip_lo = test_db->count;
          for (tmp_ctr = 0; tmp_ctr < val_db->count; tmp_ctr++)
            add_seq(test_db, val_db->list + tmp_ctr, abc);
          si.merged_skip_hi = test_db->count;
          si.merged_skip_in_test = 1;
          /* test_db now includes a val copy tail; require more total count to reach real test_limit */
          eff_test_limit += val_db->count;
        }
        val_done = TRUE;
        current_assign_frac = 0.5;
      }

      /* Housekeeping after reaching test limit */
      if ((!test_done) && (test_db->count >= eff_test_limit)) {
        if (write_test) {
          int skl = -1, skh = -1;
          if (si.merged_skip_lo >= 0 && si.merged_skip_in_test) {
            skl = si.merged_skip_lo;
            skh = si.merged_skip_hi;
          }
          save_seqs(NULL, si.test_path, si.task_id, test_db, skl, skh);
          write_test = FALSE;
        }
        /* Test reached before val: embed test into train and val. */
        if (!val_done) {
          si.train_skip_lo = train_db->count;
          for (tmp_ctr = 0; tmp_ctr < test_db->count; tmp_ctr++)
            add_seq(train_db, test_db->list + tmp_ctr, abc);
          si.train_skip_hi = train_db->count;
          si.merged_skip_lo = val_db->count;
          for (tmp_ctr = 0; tmp_ctr < test_db->count; tmp_ctr++)
            add_seq(val_db, test_db->list + tmp_ctr, abc);
          si.merged_skip_hi = val_db->count;
          si.merged_skip_in_test = 0;
          /* val_db now includes a test copy tail; require more total count to reach real val_limit */
          eff_val_limit += test_db->count;
        }
        test_done = TRUE;
      }
    }

    /* Print progress based on reporting frequency*/
    if ((seq_ctr % freq == 0) && !suppress) {

      elapsed_time = get_coarse_time() - start_time;
      print_progress(seq_ctr+1, num_records, elapsed_time);
      fflush(stderr);
    }
  }
  
  /* Print progress for one final time */
  if (!suppress) {
    elapsed_time = get_coarse_time() - start_time;
    print_progress(seq_ctr, num_records, elapsed_time);
  }
  
  fprintf(stderr, "\nEnded after going through %d sequences\n", seq_ctr);

  /* Writing train sequences to a file */
  if (write_train) {
    int skl = -1, skh = -1;
    if (si.train_skip_lo >= 0) {
      skl = si.train_skip_lo;
      skh = si.train_skip_hi;
    }
    save_seqs(NULL, si.train_path, si.task_id, train_db, skl, skh);
  }

  /* Database exhausted before test/val limits: write whatever was accumulated (same skip rules as at-limit saves). */
  if (write_test && !test_done) {
    int skl = -1, skh = -1;
    if (si.merged_skip_lo >= 0 && si.merged_skip_in_test) {
      skl = si.merged_skip_lo;
      skh = si.merged_skip_hi;
    }

    save_seqs(NULL, si.test_path, si.task_id, test_db, skl, skh);
  }
  if (write_val && !val_done) {
    int skl = -1, skh = -1;
    if (si.merged_skip_lo >= 0 && !si.merged_skip_in_test) {
      skl = si.merged_skip_lo;
      skh = si.merged_skip_hi;
    }

    save_seqs(NULL, si.val_path, si.task_id, val_db, skl, skh);
  }

  if (strcmp(si.out_path, "-") == 0)
    si.out_fp = stdout;
  else {
    len = snprintf(NULL, 0, "%s_%d%s.txt", si.out_path, si.task_id, esl_opt_GetString(go, "--shard_id")) + 1;
    ESL_ALLOC(out_file, len*sizeof(char));
    snprintf(out_file, len, "%s_%d%s.txt", si.out_path, si.task_id, esl_opt_GetString(go, "--shard_id"));
    if ((si.out_fp = fopen(out_file, "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n", out_file);
  }

  {
    int tr_emb = (si.train_skip_lo >= 0) ? (si.train_skip_hi - si.train_skip_lo) : 0;
    int mg     = (si.merged_skip_lo >= 0) ? (si.merged_skip_hi - si.merged_skip_lo) : 0;
    int logical_train = train_db->count - tr_emb;
    int logical_test  = test_db->count;
    int logical_val   = val_db->count;
    if (si.merged_skip_lo >= 0 && si.merged_skip_in_test) logical_test -= mg;
    if (si.merged_skip_lo >= 0 && !si.merged_skip_in_test) logical_val -= mg;

    fprintf(si.out_fp,
          "SledgeHMMER Stats:\n"
          "---------------------\n"
          "Train:       %d\n"
          "Test:        %d\n"
          "Val:         %d\n"
          "Discarded:   %d\n"
          "Processed:   %d\n"
          "Total:       %d\n"
          "Time:        %.2f s\n\n",
          logical_train,
          logical_test,
          logical_val,
          seq_ctr - logical_train - logical_test - logical_val,
          seq_ctr, num_records, elapsed_time);
    if (!test_done)
      fprintf(si.out_fp,
          "Note: --test_limit effective target %d was not reached (database exhausted).\n",
          eff_test_limit);
    if (!val_done)
      fprintf(si.out_fp,
          "Note: --val_limit effective target %d was not reached (database exhausted).\n",
          eff_val_limit);
  }
  
  /* Close open files */
  fclose(discard_fp);

  if (si.out_fp != stdout) {
    fclose(si.out_fp);
    free(out_file);
  }

  /* Mop the floor and run away */
  free(discard_file);
  free(results);

  esl_sq_DestroyBlock(train_db);
  esl_sq_DestroyBlock(test_db);
  esl_sq_DestroyBlock(val_db);
  esl_sq_DestroyBlock(qdb);
  esl_sq_DestroyBlock(cdb);
  esl_sq_DestroyBlock(discard_db);
  esl_sq_DestroyBlock(sqdb);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rnd);

  exit(status);

  ERROR:
    exit(status);
}