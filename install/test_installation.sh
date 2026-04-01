#!/bin/bash

set -u
set -o pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: test_installation.sh <path_to_sledge>" >&2
  echo "Example path: test_installation.sh software/sledge" >&2
  echo "Example (current directory): test_installation.sh ." >&2
  exit 2
fi
SLEDGE_DIR="$(cd "$1" && pwd)"

ROOT_DIR="${SLEDGE_DIR}/install"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
INFO_DIR="${ROOT_DIR}/info/install_${RUN_ID}"
ARTIFACT_DIR="${INFO_DIR}/artifacts"

mkdir -p "${INFO_DIR}" "${ARTIFACT_DIR}"

TEST_FASTA="${ROOT_DIR}/test_db.fasta"

SPLITTER_BIN="${SLEDGE_DIR}/bin/sledge_splitter"
PHMMER_FILTER_BIN="${SLEDGE_DIR}/bin/phmmer_filter"
SLEDGE_FILTER_BIN="${SLEDGE_DIR}/bin/sledge_filter"

TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0

# Paths for sledge_filter integration (pmbs pipeline).
# Accept MMSEQS_BIN or MMSEQS (same for BLAST_DIR_*, FASTA_DIR_*). If unset and not on PATH, values are read from test_filter.config.
# In test_filter.config, MMSEQS / BLAST_DIR / FASTA_DIR may be relative to SLEDGE_DIR (e.g. external_tools/mmseqs/bin/mmseqs).
# sledge_filter runs with cwd under install/info/.../artifacts/...; relative paths would resolve there and fail. Expand to absolute paths before export.
MMSEQS_BIN="${MMSEQS_BIN:-${MMSEQS:-}}"
BLAST_DIR_BIN="${BLAST_DIR_BIN:-${BLAST_DIR:-}}"
FASTA_DIR_BIN="${FASTA_DIR_BIN:-${FASTA_DIR:-}}"

load_filter_paths_from_config() {
  local cfg="${ROOT_DIR}/test_filter.config"
  [[ -f "${cfg}" ]] || return 1
  [[ -z "${MMSEQS_BIN}" ]] && MMSEQS_BIN="$(sed -n 's/^MMSEQS="\(.*\)"$/\1/p' "${cfg}" | head -1)"
  [[ -z "${BLAST_DIR_BIN}" ]] && BLAST_DIR_BIN="$(sed -n 's/^BLAST_DIR="\(.*\)"$/\1/p' "${cfg}" | head -1)"
  [[ -z "${FASTA_DIR_BIN}" ]] && FASTA_DIR_BIN="$(sed -n 's/^FASTA_DIR="\(.*\)"$/\1/p' "${cfg}" | head -1)"
  return 0
}

# If path is not absolute, treat it as relative to SLEDGE_DIR (matches shipped test_filter.config layout).
abs_tool_path_from_sledge() {
  local p="$1"
  [[ -z "${p}" ]] && { echo ""; return; }
  if [[ "${p}" == /* ]]; then
    echo "${p}"
  else
    echo "${SLEDGE_DIR}/${p}"
  fi
}

# Fill optional tool paths from test_filter.config, then make them cwd-independent for sledge_filter.
load_filter_paths_from_config || true
MMSEQS_BIN="$(abs_tool_path_from_sledge "${MMSEQS_BIN}")"
BLAST_DIR_BIN="$(abs_tool_path_from_sledge "${BLAST_DIR_BIN}")"
FASTA_DIR_BIN="$(abs_tool_path_from_sledge "${FASTA_DIR_BIN}")"
# run_test uses bash -lc; exported vars are visible to write_filter_config_to inside the child shell.
export SLEDGE_DIR MMSEQS_BIN BLAST_DIR_BIN FASTA_DIR_BIN PHMMER_FILTER_BIN ROOT_DIR

# Build a runnable config from test_filter.config with portable tool paths (keep in sync with test_filter.config).
write_filter_config_to() {
  local out="$1"
  awk -v pf="${PHMMER_FILTER_BIN}" -v mm="${MMSEQS_BIN}" -v bd="${BLAST_DIR_BIN}" -v fd="${FASTA_DIR_BIN}" '
    /^OUT_DIR=/ { print "OUT_DIR=\"filter_test\""; next }
    /^PHMMER_FILTER=/ { print "PHMMER_FILTER=\"" pf "\""; next }
    /^MMSEQS=/ { print "MMSEQS=\"" mm "\""; next }
    /^BLAST_DIR=/ { print "BLAST_DIR=\"" bd "\""; next }
    /^FASTA_DIR=/ { print "FASTA_DIR=\"" fd "\""; next }
    { print }
  ' "${ROOT_DIR}/test_filter.config" > "$out"
}

print_banner() {
  echo "============================================================"
  echo "Sledge installation/functionality test"
  echo "Run ID: ${RUN_ID}"
  echo "Results: ${INFO_DIR}"
  echo "Per-test logs (stdout/stderr): ${INFO_DIR}/artifacts/<NNN>_*/stdout.log and stderr.log"
  echo "============================================================"
}

check_prereqs() {
  local ok=0
  if [[ ! -d "${SLEDGE_DIR}" ]]; then
    echo "[FATAL] SLEDGE_DIR does not exist: ${SLEDGE_DIR}"
    ok=1
  fi
  if [[ ! -d "${ROOT_DIR}" ]]; then
    echo "[FATAL] Missing install directory under SLEDGE_DIR: ${ROOT_DIR}"
    ok=1
  fi
  if [[ ! -f "${TEST_FASTA}" ]]; then
    echo "[FATAL] Missing test file: ${TEST_FASTA}"
    ok=1
  fi
  if [[ ! -x "${SPLITTER_BIN}" ]]; then
    echo "[FATAL] sledge_splitter not found/executable in bin/ or src/"
    ok=1
  fi
  if [[ ! -x "${PHMMER_FILTER_BIN}" ]]; then
    echo "[FATAL] phmmer_filter not found/executable in bin/ or src/"
    ok=1
  fi
  if [[ ! -x "${SLEDGE_FILTER_BIN}" ]]; then
    echo "[FATAL] sledge_filter not found/executable in bin/"
    ok=1
  fi
  if [[ ! -f "${ROOT_DIR}/test_filter.config" ]]; then
    echo "[FATAL] Missing ${ROOT_DIR}/test_filter.config"
    ok=1
  fi
  if [[ ! -f "${ROOT_DIR}/test_db.fasta" ]]; then
    echo "[FATAL] Missing ${ROOT_DIR}/test_db.fasta"
    ok=1
  fi
  return "${ok}"
}

count_fasta() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo 0
    return 0
  fi
  local n
  n="$(grep -c '^>' "$f" 2>/dev/null || true)"
  if [[ -z "${n}" ]]; then n=0; fi
  echo "${n}"
}

is_crash_signature() {
  local f="$1"
  if [[ ! -f "$f" ]]; then return 1; fi
  if grep -Eiq '(segmentation fault|core dumped|SIGSEGV|bus error|illegal instruction|abort(ed)?)' "$f"; then
    return 0
  fi
  return 1
}

has_meaningful_error_text() {
  local out="$1"
  local err="$2"
  local txt=""
  if [[ -f "$out" ]]; then txt="${txt}$(cat "$out")"$'\n'; fi
  if [[ -f "$err" ]]; then txt="${txt}$(cat "$err")"$'\n'; fi
  if echo "$txt" | grep -Eiq '(error|failed|invalid|usage|not found|cannot|can.t|unrecognized|unexpected|unknown|must be|must set|not executable|missing|must contain|directory containing|non-negative integer|positive integer|config must)'; then
    return 0
  fi
  return 1
}

run_test() {
  local name="$1"
  local expect_fail="$2"
  local cmd="$3"
  local check_cmd="${4:-}"

  TOTAL=$((TOTAL + 1))
  local tid
  tid="$(printf "%03d" "${TOTAL}")"
  local tdir="${ARTIFACT_DIR}/${tid}_$(echo "${name}" | tr ' ' '_' | tr -cd '[:alnum:]_')"
  mkdir -p "${tdir}"

  echo ""
  echo "[TEST ${tid}] ${name}"
  echo "CMD: ${cmd}"

  bash -lc "$(declare -f write_filter_config_to); cd '${tdir}' && ${cmd}" >"${tdir}/stdout.log" 2>"${tdir}/stderr.log"
  local rc=$?

  local pass=0
  if [[ "${expect_fail}" -eq 1 ]]; then
    if [[ "${rc}" -ne 0 ]]; then
      pass=1
      if is_crash_signature "${tdir}/stdout.log" || is_crash_signature "${tdir}/stderr.log"; then
        pass=0
      fi
      if ! has_meaningful_error_text "${tdir}/stdout.log" "${tdir}/stderr.log"; then
        pass=0
      fi
    fi
  else
    [[ "${rc}" -eq 0 ]] && pass=1
  fi

  if [[ "${pass}" -eq 1 && -n "${check_cmd}" ]]; then
    bash -lc "cd '${tdir}' && $(declare -f count_fasta); ${check_cmd}" >>"${tdir}/stdout.log" 2>>"${tdir}/stderr.log"
    [[ $? -eq 0 ]] || pass=0
  fi

  if [[ "${pass}" -eq 1 ]]; then
    echo "RESULT: PASS"
    PASSED=$((PASSED + 1))
  else
    echo "RESULT: FAIL (exit=${rc})"
    echo "  If this test failed, inspect: ${tdir}/stdout.log"
    echo "  If this test failed, inspect: ${tdir}/stderr.log (intermediate configs, if any, are appended here on failure)"
    echo "  All run artifacts are under: ${INFO_DIR}"
    FAILED=$((FAILED + 1))
  fi
}

print_summary() {
  local summary="${INFO_DIR}/summary.txt"
  {
    echo "============================================================"
    echo "Run ID: ${RUN_ID}"
    echo "Total tests: ${TOTAL}"
    echo "Passed: ${PASSED}"
    echo "Failed: ${FAILED}"
    if [[ "${FAILED}" -eq 0 ]]; then
      echo "Overall: PASS"
    else
      echo "Overall: FAIL"
    fi
    echo "Skipped: ${SKIPPED}"
    echo "============================================================"
  } | tee "${summary}"
}

print_banner
check_prereqs || exit 2

# ---------------------------
# sledge_splitter tests
# ---------------------------
# Shared splitter flags (no -Z / output dir here so invalid -Z and custom-path tests can vary).
SPLIT_BASE_OPTS="--seed 42 --suppress"
# Standard runs: save default split FASTAs in each artifact test-case directory.
COMMON_SPLITTER_NO_Z="${SPLIT_BASE_OPTS} --output_dir ."
COMMON_SPLITTER="${COMMON_SPLITTER_NO_Z} -Z 500000"

run_test \
  "splitter single cpu single init_chunk" \
  0 \
  "'${SPLITTER_BIN}' ${COMMON_SPLITTER} --cpu 1 --dbblock 100 --init_chunk 1 --test_limit 20 --val_limit 10 -o stats '${TEST_FASTA}'" \
  "DB_N=\$(count_fasta '${TEST_FASTA}'); \
   TRAIN_N=\$(count_fasta train_0.fasta); \
   TEST_N=\$(count_fasta test_0.fasta); \
   VAL_N=\$(count_fasta val_0.fasta); \
   PROC_N=\$(awk '/^Processed:/ {print \$2}' stats_0.txt); \
   test \"\${TEST_N}\" -ge 20 && test \"\${VAL_N}\" -ge 10 && \
   test -n \"\${PROC_N}\" && test \"\${PROC_N}\" -le \"\${DB_N}\" && \
   test \$((TRAIN_N + TEST_N + VAL_N)) -eq \"\${PROC_N}\""

run_test \
  "splitter multi cpu multi init_chunk" \
  0 \
  "'${SPLITTER_BIN}' ${COMMON_SPLITTER} --cpu 8 --dbblock 100 --init_chunk 8 --test_limit 20 --val_limit 10 -o stats '${TEST_FASTA}'" \
  "DB_N=\$(count_fasta '${TEST_FASTA}'); \
   TRAIN_N=\$(count_fasta train_0.fasta); \
   TEST_N=\$(count_fasta test_0.fasta); \
   VAL_N=\$(count_fasta val_0.fasta); \
   PROC_N=\$(awk '/^Processed:/ {print \$2}' stats_0.txt); \
   test \"\${TEST_N}\" -ge 20 && test \"\${VAL_N}\" -ge 10 && \
   test -n \"\${PROC_N}\" && test \"\${PROC_N}\" -le \"\${DB_N}\" && \
   test \$((TRAIN_N + TEST_N + VAL_N)) -eq \"\${PROC_N}\""

run_test \
  "splitter test limit reached before val limit" \
  0 \
  "'${SPLITTER_BIN}' ${COMMON_SPLITTER} --cpu 1 --dbblock 100 --init_chunk 4 --test_limit 10 --val_limit 20 -o stats '${TEST_FASTA}'" \
  "DB_N=\$(count_fasta '${TEST_FASTA}'); \
   TRAIN_N=\$(count_fasta train_0.fasta); \
   TEST_N=\$(count_fasta test_0.fasta); \
   VAL_N=\$(count_fasta val_0.fasta); \
   PROC_N=\$(awk '/^Processed:/ {print \$2}' stats_0.txt); \
   test \"\${TEST_N}\" -ge 10 && test \"\${VAL_N}\" -ge 20 && \
   test -n \"\${PROC_N}\" && test \"\${PROC_N}\" -le \"\${DB_N}\" && \
   test \$((TRAIN_N + TEST_N + VAL_N)) -eq \"\${PROC_N}\""

run_test \
  "splitter disable val mode" \
  0 \
  "'${SPLITTER_BIN}' ${COMMON_SPLITTER} --disable_val --cpu 1 --dbblock 100 --init_chunk 4 --test_limit 25 -o stats '${TEST_FASTA}'" \
  "DB_N=\$(count_fasta '${TEST_FASTA}'); \
   TRAIN_N=\$(count_fasta train_0.fasta); \
   TEST_N=\$(count_fasta test_0.fasta); \
   PROC_N=\$(awk '/^Processed:/ {print \$2}' stats_0.txt); \
   test ! -f val_0.fasta && test \"\${TEST_N}\" -ge 25 && \
   test -n \"\${PROC_N}\" && test \"\${PROC_N}\" -le \"\${DB_N}\" && \
   test \$((TRAIN_N + TEST_N)) -eq \"\${PROC_N}\""

run_test \
  "splitter invalid input path error handling" \
  1 \
  "'${SPLITTER_BIN}' ${COMMON_SPLITTER} --cpu 1 --dbblock 100 --init_chunk 1 --test_limit 10 --val_limit 5 -o stats '${ROOT_DIR}/missing_input.fasta'"

run_test \
  "splitter invalid db size option handling" \
  1 \
  "'${SPLITTER_BIN}' ${COMMON_SPLITTER_NO_Z} -Z -3 --cpu 1 --dbblock 100 --init_chunk 1 --test_limit 10 --val_limit 5 -o stats '${TEST_FASTA}'"

run_test \
  "splitter custom train test val discard path prefixes" \
  0 \
  "mkdir -p custom_out && '${SPLITTER_BIN}' ${SPLIT_BASE_OPTS} -Z 500000 --cpu 1 --dbblock 100 --init_chunk 1 --test_limit 20 --val_limit 10 --train_path custom_out/mytrain --test_path custom_out/mytest --val_path custom_out/myval --discard_path custom_out/mydiscard -o stats '${TEST_FASTA}'" \
  "DB_N=\$(count_fasta '${TEST_FASTA}'); \
   TRAIN_N=\$(count_fasta custom_out/mytrain_0.fasta); \
   TEST_N=\$(count_fasta custom_out/mytest_0.fasta); \
   VAL_N=\$(count_fasta custom_out/myval_0.fasta); \
   PROC_N=\$(awk '/^Processed:/ {print \$2}' stats_0.txt); \
   test -f custom_out/mydiscard_0.fasta && \
   test \"\${TEST_N}\" -ge 20 && test \"\${VAL_N}\" -ge 10 && \
   test -n \"\${PROC_N}\" && test \"\${PROC_N}\" -le \"\${DB_N}\" && \
   test \$((TRAIN_N + TEST_N + VAL_N)) -eq \"\${PROC_N}\""

# ---------------------------
# sledge_filter tests
# ---------------------------
# Base without -Z / qblock / tblock so invalid -Z and undersized block tests set each once.
COMMON_FILTER_BASE="--seed 42 --suppress"
COMMON_FILTER_NO_Z="${COMMON_FILTER_BASE} --qblock 100 --tblock 100"
COMMON_FILTER="${COMMON_FILTER_NO_Z} -Z 500000"

run_test \
  "filter format 0 single cpu single query per core" \
  0 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --cpu 1 --qsize 1 --format 0 -o filt_fmt0 '${TEST_FASTA}' '${TEST_FASTA}'" \
  "test -f filt_fmt0_0.txt"

run_test \
  "filter format 1 multi cpu single query per core" \
  0 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --cpu 4 --qsize 1 --format 1 -o filt_fmt1 '${TEST_FASTA}' '${TEST_FASTA}'" \
  "test -s filt_fmt1_0.txt"

run_test \
  "filter format 2 multi cpu multi query per core" \
  0 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --cpu 1 --qsize 5 --format 2 -o filt_fmt2 '${TEST_FASTA}' '${TEST_FASTA}'" \
  "test -f filt_fmt2_0.txt"

run_test \
  "filter format 3 reject query ids" \
  0 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --cpu 1 --qsize 4 --format 3 -o filt_fmt3 '${TEST_FASTA}' '${TEST_FASTA}'" \
  "test -f filt_fmt3_0.txt"

run_test \
  "filter format 4 reject target ids all_hits mode" \
  0 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --all_hits --cpu 1 --qsize 3 --format 4 -o filt_fmt4 '${TEST_FASTA}' '${TEST_FASTA}'" \
  "test -f filt_fmt4_0.txt"

run_test \
  "filter invalid input path error handling" \
  1 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --cpu 1 --qsize 1 --format 0 -o filt_bad '${ROOT_DIR}/missing_q.fasta' '${TEST_FASTA}'"

run_test \
  "filter invalid declared format handling" \
  1 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER} --cpu 1 --qsize 1 --format 0 --qformat definitely_not_a_format -o filt_badfmt '${TEST_FASTA}' '${TEST_FASTA}'"

run_test \
  "filter invalid db size option handling" \
  1 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER_NO_Z} -Z 0 --cpu 1 --qsize 1 --format 0 -o filt_badz '${TEST_FASTA}' '${TEST_FASTA}'"

run_test \
  "filter undersized qblock graceful failure" \
  1 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER_BASE} -Z 500000 --cpu 1 --qsize 1 --format 0 --qblock 5 --tblock 100 -o filt_small_qblock '${TEST_FASTA}' '${TEST_FASTA}'"

run_test \
  "filter undersized tblock graceful failure" \
  1 \
  "'${PHMMER_FILTER_BIN}' ${COMMON_FILTER_BASE} -Z 500000 --cpu 1 --qsize 1 --format 0 --qblock 100 --tblock 5 -o filt_small_tblock '${TEST_FASTA}' '${TEST_FASTA}'"

# ---------------------------
# sledge_filter (pmbs pipeline) integration
# ---------------------------
run_test \
  "sledge_filter pipeline pmbs install integration" \
  0 \
  "write_filter_config_to ft.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}' --out-suffix check" \
  "test -d filter_test/phmmer_check && test -d filter_test/mmseqs_check && test -d filter_test/blast_check && test -d filter_test/sw_check"

run_test \
  "sledge_filter pipeline spm order" \
  0 \
  "write_filter_config_to ft.config && '${SLEDGE_FILTER_BIN}' --order spm --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}' --out-suffix check" \
  "test -d filter_test/phmmer_check && test -d filter_test/mmseqs_check && test -d filter_test/sw_check"

run_test \
  "sledge_filter missing required args" \
  1 \
  "'${SLEDGE_FILTER_BIN}'"

run_test \
  "sledge_filter missing config file" \
  1 \
  "'${SLEDGE_FILTER_BIN}' --order pmbs --config ./does_not_exist_$$.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter empty config missing OUT_DIR" \
  1 \
  ": > empty.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./empty.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter unknown order character" \
  1 \
  "'${SLEDGE_FILTER_BIN}' --order pmxs --config ./irrelevant.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter invalid REMOVE_TARGET in config" \
  1 \
  "write_filter_config_to ft.config && sed 's/^REMOVE_TARGET=.*/REMOVE_TARGET=\"bad\"/' ft.config > ft2.config && mv ft2.config ft.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter invalid SLEDGE_CORES in config" \
  1 \
  "write_filter_config_to ft.config && sed 's/^SLEDGE_CORES=.*/SLEDGE_CORES=\"not_a_number\"/' ft.config > ft2.config && mv ft2.config ft.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter invalid SLEDGE_FILTER path in config" \
  1 \
  "write_filter_config_to ft.config && sed 's|^PHMMER_FILTER=.*|PHMMER_FILTER=\"/nonexistent/phmmer_filter\"|' ft.config > ft2.config && mv ft2.config ft.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter invalid MMSEQS path in config" \
  1 \
  "write_filter_config_to ft.config && sed 's|^MMSEQS=.*|MMSEQS=\"/nonexistent/mmseqs\"|' ft.config > ft2.config && mv ft2.config ft.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

run_test \
  "sledge_filter invalid Z_SIZE in config" \
  1 \
  "write_filter_config_to ft.config && sed 's/^Z_SIZE=.*/Z_SIZE=\"0\"/' ft.config > ft2.config && mv ft2.config ft.config && '${SLEDGE_FILTER_BIN}' --order pmbs --config ./ft.config --fixed-file '${ROOT_DIR}/test_fixed.fasta' --db-file '${TEST_FASTA}'"

print_summary

if [[ "${FAILED}" -eq 0 ]]; then
  exit 0
fi
exit 1
