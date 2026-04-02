/* sledge_dev_og.h — public API for sledge_dev_og.c (esl_threads worker model).
 *
 * Derived from sledge_dev.h: same structs and entry points as the current
 * splitter/phmmer drivers expect, except:
 *   - store_results matches the OG implementation (takes ESL_THREADS *).
 *   - sledge_worker_pool_destroy is omitted (OG does not use a persistent pool).
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_threads.h"

#ifdef HMMER_THREADS
#include <unistd.h>
#endif

#include "hmmer.h"

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 *
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char  *qfile;        /* query sequence file                                   */
  char  *dbfile;       /* database file                                         */

  int    do_mpi;       /* TRUE if we're doing MPI parallelization               */
  int    nproc;        /* how many MPI processes, total                         */
  int    my_rank;      /* who am I, in 0..nproc-1                               */

  char  *firstseq_key; /* name of the first sequence in the restricted db range */
  int    n_targetseq;  /* number of sequences in the restricted range           */
};

typedef struct {
  char   decision;     /* Store accept or reject       */
  float *pids;         /* List of failed hit PIDS      */
  char **targets;      /* Non-owning target IDs        */
  char  *source;       /* Non-owning query ID          */
  int   *lengths;      /* Ali length of failed hits    */
  float *evals;        /* E-value of failed hits       */
  int    n_targets;    /* Number of failed hits        */
  int    capacity;     /* Capacity of faield hits      */
  float  max_pid;      /* Max PID of accepted sequence */
} RESULT_INFO;

typedef struct {
  int    qsize;        /* Number of queries per thread            */
  int    cores;        /* Number of cpu cores to use              */
  float  pid_low;      /* Lower limit for PID threshold           */
  float  pid_high;     /* Upper limit for PID threshold           */
  float  E;            /* Inclusion E-Value threshold             */
  int    db_size;      /* Size of target database                 */
  bool   early_stop;   /* Flag for early stopping                 */
  int    format;       /* Output format for filter                */
  bool   train_only;   /* Flag for train_only mode                */
  bool   test_only;    /* Flag for test_only mode                 */
  bool   low_only;     /* Flag to stop after test/val limit       */
  float  train_frac;   /* Max train_fraction                      */
  int    init_chunk;   /* Initial query chunk size                */
  int    test_limit;   /* Minimum required test sequences         */
  int    val_limit;    /* Minimum required val sequences          */
  int    task_id;      /* Task ID or shard number for process     */
  char  *train_path;   /* File path for train db                  */
  char  *test_path;    /* File path for test db                   */
  char  *val_path;     /* File path for val db                    */
  char  *discard_path; /* File path for discard db                */
  char  *out_path;     /* File path for output                    */
  FILE  *out_fp;       /* File ptr for output                     */
  int    offset;       /* Track current buffer size               */
  int    train_skip_lo; /* if >= 0, omit [lo, hi) from train FASTA (embedded copy) */
  int    train_skip_hi;
  int    merged_skip_lo; /* embedded copy in test or val FASTA; see merged_skip_in_test */
  int    merged_skip_hi;
  int    merged_skip_in_test; /* 1=merged range is in test_db, 0=merged range is in val_db */
} SLEDGE_INFO;

typedef struct {
  P7_BUILDER       *bld;         /* p7 builder object             */
  ESL_GETOPTS      *go;          /* Commandline options           */
  RESULT_INFO      *result;      /* Store accept/reject results   */
  int               result_cap;  /* Allocated RESULT_INFO slots   */
  ESL_SQ_BLOCK     *sqdb;        /* Pointer to sequence database  */
  ESL_SQ_BLOCK     *qdb;         /* Pointer to query database     */
  SLEDGE_INFO      *si;          /* Pointer to sledgehmmer info   */
  int               start;       /* Query start index             */
  int               num_queries; /* Number of queries to handle   */
} WORKER_INFO;

/* License information */
#define SHMMER_VERSION "1.0"
#define SHMMER_DATE "May 2025"
#define SHMMER_COPYRIGHT "Copyright (C) 2025 Howard Hughes Medical Institute."
#define SHMMER_LICENSE "Freely distributed under the BSD open source license."
#define SHMMER_URL "http://hmmer.org/"

/* p7 pipeline options */
#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#if defined (HMMER_THREADS) && defined (HMMER_MPI)
  #define CPUOPTS     "--mpi"
  #define MPIOPTS     "--cpu"
#else
  #define CPUOPTS     NULL
  #define MPIOPTS     NULL
#endif

/* Sledgehmmer defaults */
#define ACCEPT '1'
#define REJECT '0'
#define INITIAL_HITS_CAPACITY 20
#define BUFFER_MAX 65536

/* Sledgehmmer output formats */
#define R_STRING  0 /* ACCEPT/REJECT string     */
#define R_ALL     1 /* Full info of reject hits */
#define R_ACCEPT  2 /* IDs of ACCEPT sequences  */
#define R_REJECTQ 3 /* IDs of REJECT queries    */
#define R_REJECTT 4 /* IDs of REJECT targets    */

extern void sledge_banner(FILE *fp, const char *progname, char *banner);

extern int read_cust(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock);

extern int destroy_info(WORKER_INFO *info, int num_cores);

static inline int append_target(RESULT_INFO *result, const char *target_id, float pid, float eval, int length);

extern int add_seq(ESL_SQ_BLOCK *db, ESL_SQ *seq, const ESL_ALPHABET *abc);

extern void add_to_buffer(char *buffer, char *line, size_t line_length, SLEDGE_INFO *si);

extern int store_results(ESL_THREADS *threadObj, SLEDGE_INFO *si, char *results);

extern double get_coarse_time(void);

extern void format_time(double seconds, char *buffer, size_t buffer_size);

extern void print_progress(int current, int total, double elapsed);

extern int save_seqs(FILE *fp, char *fpath, int fid, ESL_SQ_BLOCK *db, int skip_lo, int skip_hi);

extern int check_triple(ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *db_1, ESL_SQ_BLOCK *db_2,
    ESL_SQ_BLOCK *db_3, FILE *discard_fp, SLEDGE_INFO *si, ESL_ALPHABET *abc,
    ESL_GETOPTS *go);

extern int process_two_group_chunk(ESL_GETOPTS *go, ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *db_a,
    ESL_SQ_BLOCK *db_b, ESL_SQ_BLOCK *check_a, ESL_SQ_BLOCK *check_b,
    ESL_SQ_BLOCK *cdb, ESL_SQ_BLOCK *discard_db, FILE *discard_fp,
    ESL_ALPHABET *abc, SLEDGE_INFO *si, bool assign_a, char *results);

extern int assign_master(ESL_GETOPTS *go, ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *sqdb, SLEDGE_INFO *si, char *results);

extern int p7_Pipeline_cust(P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq,
  const ESL_SQ *qsq, P7_TOPHITS *hitlist, RESULT_INFO *result, SLEDGE_INFO *si);

extern int p7_domaindef_cust(const ESL_SQ *sq, const ESL_SQ *qsq, P7_OPROFILE *om,
  P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck,
  P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
  P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr,
  RESULT_INFO *result, float *nullsc, SLEDGE_INFO *si);

static int is_multidomain_region(P7_DOMAINDEF *ddef, int i, int j);

static int region_trace_ensemble  (P7_DOMAINDEF *ddef, const P7_OPROFILE *om, const ESL_DSQ *dsq, int ireg, int jreg,
  const P7_OMX *fwd, P7_OMX *wrk, int *ret_nc);

static int rescore_domain_cust(P7_DOMAINDEF *ddef, P7_OPROFILE *om, const ESL_SQ *sq, const ESL_SQ *qsq,
  P7_OMX *ox1, P7_OMX *ox2, int i, int j, int null2_is_done, P7_BG *bg, int long_target,
  P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr, RESULT_INFO *result,
  float *nullsc, SLEDGE_INFO *si);

static int find_pid(P7_TRACE *tr, const ESL_SQ *sq, const ESL_SQ *qsq, float *pid, char *decision, SLEDGE_INFO *si);

static void pipeline_thread(void *arg);