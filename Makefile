SHELL := /bin/sh

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD := $(ROOT)/build
BIN := $(ROOT)/bin

UNAME_S := $(shell uname -s 2>/dev/null || echo Unknown)
OS ?= $(UNAME_S)
EXE :=
ifeq ($(OS),Windows_NT)
  EXE := .exe
endif

CC ?= gcc
AR ?= ar
RANLIB ?= ranlib

CFLAGS := -O3 -pthread
CPPFLAGS := -DHAVE_CONFIG_H
INCLUDES := -I$(ROOT)/easel -I$(ROOT)/src -I$(ROOT)/src/impl_sse
LDFLAGS :=
LDLIBS := -lpthread -lm

SLEDGE_SRCS := \
	src/sledge_splitter.c \
	src/phmmer_filter.c \
	src/sledge_dev.c

HMMER_SRCS := \
	src/errors.c src/logsum.c src/p7_alidisplay.c src/p7_bg.c src/p7_builder.c src/p7_domaindef.c \
	src/p7_hmm.c src/p7_pipeline.c src/p7_prior.c src/p7_profile.c src/p7_spensemble.c src/p7_tophits.c \
	src/p7_trace.c src/p7_scoredata.c src/fm_general.c src/fm_sse.c src/fm_ssv.c src/build.c src/evalues.c \
	src/eweight.c src/hmmer.c src/modelconfig.c src/modelstats.c src/seqmodel.c src/tracealign.c src/p7_gmx.c \
	src/p7_hit.c src/p7_hmmwindow.c src/fm_alphabet.c src/emit.c src/generic_decoding.c src/generic_fwdback.c \
	src/generic_optacc.c src/p7_domain.c src/impl_sse/decoding.c src/impl_sse/fwdback.c src/impl_sse/io.c \
	src/impl_sse/msvfilter.c src/impl_sse/null2.c src/impl_sse/optacc.c src/impl_sse/stotrace.c \
	src/impl_sse/vitfilter.c src/impl_sse/p7_omx.c src/impl_sse/p7_oprofile.c src/impl_sse/ssvfilter.c

EASEL_SRCS := \
	easel/easel.c easel/esl_alphabet.c easel/esl_cluster.c easel/esl_dirichlet.c easel/esl_dmatrix.c easel/esl_exponential.c \
	easel/esl_fileparser.c easel/esl_getopts.c easel/esl_gumbel.c easel/esl_hmm.c easel/esl_keyhash.c easel/esl_mem.c \
	easel/esl_minimizer.c easel/esl_mixdchlet.c easel/esl_msa.c easel/esl_msacluster.c easel/esl_msaweight.c \
	easel/esl_quicksort.c easel/esl_random.c easel/esl_rand64.c easel/esl_randomseq.c easel/esl_rootfinder.c \
	easel/esl_scorematrix.c easel/esl_sq.c easel/esl_sqio.c easel/esl_sqio_ascii.c easel/esl_sqio_ncbi.c easel/esl_ssi.c \
	easel/esl_stats.c easel/esl_stopwatch.c easel/esl_threads.c easel/esl_tree.c easel/esl_vectorops.c easel/esl_wuss.c \
	easel/esl_sse.c easel/esl_arr2.c easel/esl_arr3.c easel/esl_bitfield.c easel/esl_composition.c easel/esl_distance.c \
	easel/esl_graph.c easel/esl_matrixops.c easel/esl_msafile.c easel/esl_msafile_a2m.c easel/esl_msafile_afa.c \
	easel/esl_msafile_clustal.c easel/esl_msafile_phylip.c easel/esl_msafile_psiblast.c easel/esl_msafile_selex.c \
	easel/esl_msafile_stockholm.c easel/esl_ratematrix.c easel/esl_stack.c easel/esl_buffer.c

SLEDGE_OBJS := $(patsubst %.c,$(BUILD)/%.o,$(SLEDGE_SRCS))
HMMER_OBJS := $(patsubst %.c,$(BUILD)/%.o,$(HMMER_SRCS))
EASEL_OBJS := $(patsubst %.c,$(BUILD)/%.o,$(EASEL_SRCS))

.PHONY: all libs clean install-external

all: $(BIN)/sledge_splitter$(EXE) $(BIN)/phmmer_filter$(EXE) $(BIN)/sledge_filter

libs: $(BUILD)/libhmmer_min.a $(BUILD)/libeasel_min.a

$(BUILD)/libhmmer_min.a: $(HMMER_OBJS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(BUILD)/libeasel_min.a: $(EASEL_OBJS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(BIN)/sledge_splitter$(EXE): $(BUILD)/src/sledge_splitter.o $(BUILD)/src/sledge_dev.o $(BUILD)/libhmmer_min.a $(BUILD)/libeasel_min.a
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BIN)/phmmer_filter$(EXE): $(BUILD)/src/phmmer_filter.o $(BUILD)/src/sledge_dev.o $(BUILD)/libhmmer_min.a $(BUILD)/libeasel_min.a
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BIN)/sledge_filter: $(ROOT)/src/sledge_filter.sh
	@mkdir -p $(dir $@)
	cp $< $@
	chmod +x $@

$(BUILD)/%.o: $(ROOT)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

install-external:
	"$(ROOT)/install/install_external.sh" "$(ROOT)"

clean:
	rm -rf $(BUILD) $(BIN)
