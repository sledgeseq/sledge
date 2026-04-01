#!/bin/bash
#
# filter_revised.sh — Run pHMMER, MMseqs2, BLAST, and/or SW filters in a chosen order.
# Required: --order, --config, --fixed-file, --db-file. Other parameters come from the config file.
#
# Usage:
#   --order STRING           Order of tools: p=phmmer, m=mmseqs, b=blast, s=sw (e.g. pmb, spm)
#   --config FILE            Path to config file (required)
#   --fixed-file FILE        Fixed/reference FASTA (e.g. test set)
#   --db-file FILE           DB FASTA (shard / sequence database to be filtered against fixed file)
#   --out-suffix SUFFIX      Optional. Suffix for tool output dirs (phmmer_<suffix>, mmseqs_<suffix>, etc.).
#                            Default: TASK_ID from config. Use to separate outputs per run (e.g. 0_vs_0_0).
#
# Example: ./filter_revised.sh --order pmb --config filter_revised.config --fixed-file ../data/shards/shard_0.fasta --db-file ../data/distant_test.fasta
#
# Set OUT_DIR, KEEP_INTERMEDIATES, etc. in the config file.

# #SBATCH -J filter_all
# #SBATCH -p eddy
# #SBATCH -c 96
# #SBATCH -t "15-00:00"
# #SBATCH --mem-per-cpu="1GB"
# #SBATCH --array=0
# #SBATCH -o logs/filter_revised_%a.out
# #SBATCH -e logs/filter_revised_%a.err
# #SBATCH --mem-per-cpu=1GB # Based on MaxRSS value, 20GB is total limit

set -euo pipefail

# --- Argument parsing: required --order, --config, --fixed-file, --db-file; optional --out-suffix ---
ORDER=""
CONFIG_FILE=""
FIXED_FILE=""
DB_FILE=""
OUT_SUFFIX_ARG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --order)
      ORDER="$2"
      shift 2
      ;;
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --fixed-file)
      FIXED_FILE="$2"
      shift 2
      ;;
    --db-file)
      DB_FILE="$2"
      shift 2
      ;;
    --out-suffix)
      OUT_SUFFIX_ARG="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$ORDER" || -z "$CONFIG_FILE" || -z "$FIXED_FILE" || -z "$DB_FILE" ]]; then
  echo "Usage: $0 --order STRING --config FILE --fixed-file FILE --db-file FILE"
  echo "  Order: p=phmmer, m=mmseqs, b=blast, s=sw (e.g. pmb, spm)"
  exit 1
fi

# --- Validate --order before loading config (clear errors without needing valid tool paths) ---
for (( _oi=0; _oi<${#ORDER}; _oi++ )); do
  case "${ORDER:${_oi}:1}" in
    p|P|m|M|b|B|s|S) ;;
    *)
      echo "Unknown tool in --order: ${ORDER:${_oi}:1}" >&2
      exit 1
      ;;
  esac
done

# --- Load config ---
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi
source "$CONFIG_FILE"

# Required from config
for var in OUT_DIR; do
  if [[ -z "${!var:-}" ]]; then
    echo "Config must set $var" >&2
    exit 1
  fi
done

# Optional defaults from config if not set (KEEP_INTERMEDIATES is in config only)
KEEP_INTERMEDIATES="${KEEP_INTERMEDIATES:-}"

# TASK_ID resolution:
# - If config sets TASK_ID=SLURM, use $SLURM_ARRAY_TASK_ID when available (else 0).
# - Otherwise, use the numeric TASK_ID from config, defaulting to 0 when unset.
if [[ "${TASK_ID:-}" == "SLURM" ]]; then
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    TASK_ID="${SLURM_ARRAY_TASK_ID}"
  else
    TASK_ID="0"
  fi
else
  TASK_ID="${TASK_ID:-0}"
fi

# Tool output dirs use OUT_SUFFIX (default TASK_ID); TASK_ID stays integer for internal use (e.g. sledge --task_id)
OUT_SUFFIX="${OUT_SUFFIX_ARG:-${OUT_SUFFIX:-$TASK_ID}}"
REMOVE_TARGET="${REMOVE_TARGET:-db}"
E_VALUE="${E_VALUE:-0.01}"
Z_SIZE="${Z_SIZE:-81514348}"
SLEDGE_CORES="${SLEDGE_CORES:-96}"
MMSEQS_CORES="${MMSEQS_CORES:-96}"
BLAST_CORES="${BLAST_CORES:-96}"
SW_CORES="${SW_CORES:-96}"
BLAST_DBSIZE="${BLAST_DBSIZE:-32596740121}"
BLAST_MAX_TARGET_SEQS="${BLAST_MAX_TARGET_SEQS:-12000000}"
MMSEQS_MAX_SEQS="${MMSEQS_MAX_SEQS:-12000000}"

if [[ "$REMOVE_TARGET" != "db" && "$REMOVE_TARGET" != "fixed" ]]; then
  echo "REMOVE_TARGET must be 'db' or 'fixed'" >&2
  exit 1
fi

is_uint() { [[ "${1:-}" =~ ^[0-9]+$ ]]; }
is_positive_uint() { [[ "${1:-}" =~ ^[1-9][0-9]*$ ]]; }

if ! is_positive_uint "$Z_SIZE"; then
  echo "Z_SIZE must be a positive integer (got: ${Z_SIZE:-})" >&2
  exit 1
fi
for _cname in SLEDGE_CORES MMSEQS_CORES BLAST_CORES SW_CORES; do
  _cv="${!_cname}"
  if ! is_uint "${_cv}"; then
    echo "${_cname} must be a non-negative integer (got: ${_cv:-})" >&2
    exit 1
  fi
done

# Only require binaries for tools listed in --order (e.g. pms skips BLAST if b is absent).
NEED_PHMMER=0
NEED_MMSEQS=0
NEED_BLAST=0
NEED_SW=0
for (( _ni=0; _ni<${#ORDER}; _ni++ )); do
  case "${ORDER:${_ni}:1}" in
    p|P) NEED_PHMMER=1 ;;
    m|M) NEED_MMSEQS=1 ;;
    b|B) NEED_BLAST=1 ;;
    s|S) NEED_SW=1 ;;
  esac
done

if [[ "${NEED_PHMMER}" -eq 1 ]]; then
  if [[ ! -x "${PHMMER_FILTER}" ]]; then
    echo "PHMMER_FILTER is not executable or missing: ${PHMMER_FILTER}" >&2
    exit 1
  fi
fi
if [[ "${NEED_MMSEQS}" -eq 1 ]]; then
  if [[ ! -x "${MMSEQS}" ]]; then
    echo "MMSEQS is not executable or missing: ${MMSEQS}" >&2
    exit 1
  fi
  _bc_exe="$(command -v bc 2>/dev/null || true)"
  [[ -z "${_bc_exe}" && -x /usr/bin/bc ]] && _bc_exe=/usr/bin/bc
  [[ -z "${_bc_exe}" && -x /bin/bc ]] && _bc_exe=/bin/bc
  if [[ -z "${_bc_exe}" || ! -x "${_bc_exe}" ]]; then
    echo "MMseqs step requires bc (not found on PATH or /usr/bin/bc, /bin/bc)" >&2
    exit 1
  fi
fi
if [[ "${NEED_BLAST}" -eq 1 ]]; then
  if [[ ! -d "${BLAST_DIR}" ]] || [[ ! -x "${BLAST_DIR}/makeblastdb" ]] || [[ ! -x "${BLAST_DIR}/blastp" ]]; then
    echo "BLAST_DIR must be a directory containing makeblastdb and blastp: ${BLAST_DIR}" >&2
    exit 1
  fi
fi
if [[ "${NEED_SW}" -eq 1 ]]; then
  if [[ ! -d "${FASTA_DIR}" ]] || [[ ! -x "${FASTA_DIR}/ssearch36" ]]; then
    echo "FASTA_DIR must be a directory containing ssearch36: ${FASTA_DIR}" >&2
    exit 1
  fi
fi

mkdir -p "$OUT_DIR"

# --- Reusable: remove sequences by ID list (pure bash/awk, no Python) ---
# remove_seqs_bash <input_fasta> <omit_ids_file> <output_fasta> [--long]
remove_seqs_bash() {
  local input_fasta="$1"
  local id_file="$2"
  local output_fasta="$3"
  local long_opt="${4:-}"
  if [[ ! -f "$input_fasta" ]]; then
    echo "remove_seqs_bash: input not found: $input_fasta" >&2
    return 1
  fi
  if [[ ! -f "$id_file" ]]; then
    echo "remove_seqs_bash: id file not found: $id_file" >&2
    return 1
  fi
  # If no IDs are listed to omit, keep the FASTA unchanged.
  if ! awk 'NF { found=1; exit } END { exit(found ? 0 : 1) }' "$id_file"; then
    cp "$input_fasta" "$output_fasta"
    return $?
  fi
  local use_long=0
  [[ "$long_opt" == "--long" ]] && use_long=1
  awk -v id_file="$id_file" -v use_long="$use_long" '
    BEGIN {
      while ((getline line < id_file) > 0) {
        sub(/[ \t\r]+$/, "", line)
        if (line != "") omit[line] = 1
      }
      close(id_file)
      rec = ""
      id = ""
    }
    /^>/ {
      if (rec != "" && !(id in omit)) printf "%s\n", rec
      rec = $0
      id = substr($0, 2)
      if (use_long && length(id) >= 3 && substr(id, 3, 1) == "|") {
        n = split(id, a, "|")
        id = (n >= 2 ? a[2] : id)
      } else {
        sub(/[ \t].*$/, "", id)
      }
      next
    }
    { rec = rec (rec == "" ? "" : "\n") $0 }
    END {
      if (rec != "" && !(id in omit)) printf "%s\n", rec
    }
  ' "$input_fasta" > "$output_fasta"
}

# --- Exit if a tool's output FASTA has no sequences ---
check_fasta_nonempty() {
  local fasta="$1"
  local tool_label="$2"
  local n
  n=$(grep -c ">" "$fasta" 2>/dev/null) || true
  n=${n:-0}
  if [[ "$n" -eq 0 ]]; then
    echo "No sequence survived the filtering (after $tool_label)."
    exit 0
  fi
}

# --- PHMMER ---
run_phmmer() {
  local start=$SECONDS
  local out_p="${OUT_DIR}/phmmer_${OUT_SUFFIX}"
  mkdir -p "$out_p"
  local omit_ids="${out_p}/phmmer_hits_${TASK_ID}.txt"
  : > "$omit_ids"
  local phmmer_out="${out_p}/phmmer_${TASK_ID}.fasta"

  append_omit_ids() {
    local filt_file="$1"
    local col="${2:-1}"
    awk -v c="$col" 'NR>2 {print $c}' "$filt_file" >> "$omit_ids"
  }

  local n_db n_fixed
  n_db=$(grep -c ">" "$DB_FILE" || true)
  n_fixed=$(grep -c ">" "$FIXED_FILE" || true)
  # First call: query=DB, target=FIXED → qblock=n_db, tblock=n_fixed
  $PHMMER_FILTER --cpu "$SLEDGE_CORES" --qsize 100 --qblock "$n_db" --tblock "$n_fixed" --phigh 0 --plow 0 \
    -E "$E_VALUE" -Z "$Z_SIZE" --task_id "$TASK_ID" -o "${out_p}/phmmer_hits_db_fi" \
    "$DB_FILE" "$FIXED_FILE"
  # Second call: query=FIXED, target=DB → qblock=n_fixed, tblock=n_db
  $PHMMER_FILTER --cpu "$SLEDGE_CORES" --qsize 100 --qblock "$n_fixed" --tblock "$n_db" --phigh 0 --plow 0 \
    -E "$E_VALUE" -Z "$Z_SIZE" --task_id "$TASK_ID" -o "${out_p}/phmmer_hits_fi_db" \
    "$FIXED_FILE" "$DB_FILE"


  if [[ "$REMOVE_TARGET" == "fixed" ]]; then
    append_omit_ids "${out_p}/phmmer_hits_db_fi_${TASK_ID}.txt" 2
    append_omit_ids "${out_p}/phmmer_hits_fi_db_${TASK_ID}.txt" 1
    remove_seqs_bash "$FIXED_FILE" "$omit_ids" "${out_p}/fixed_phmmer_${TASK_ID}.fasta"
    FIXED_FILE="${out_p}/fixed_phmmer_${TASK_ID}.fasta"
  else
    append_omit_ids "${out_p}/phmmer_hits_db_fi_${TASK_ID}.txt" 1
    append_omit_ids "${out_p}/phmmer_hits_fi_db_${TASK_ID}.txt" 2
    remove_seqs_bash "$DB_FILE" "$omit_ids" "$phmmer_out"
    DB_FILE="$phmmer_out"
  fi

  if [[ -z "$KEEP_INTERMEDIATES" ]]; then
    rm -f "${out_p}/results_phmmer_db_fi_${TASK_ID}.txt" "${out_p}/results_phmmer_fi_db_${TASK_ID}.txt"
  fi
  if [[ "$REMOVE_TARGET" == "db" ]]; then
    check_fasta_nonempty "$phmmer_out" "pHMMER"
  else
    check_fasta_nonempty "$FIXED_FILE" "pHMMER"
  fi
  echo "PHMMER took $((SECONDS - start)) seconds"
}

# --- MMSEQS2 (iterative) ---
run_mmseqs() {
  local start=$SECONDS
  local out_mm="${OUT_DIR}/mmseqs_${OUT_SUFFIX}"
  mkdir -p "$out_mm"
  local iteration=0
  local mm_db_file="$DB_FILE"
  local long_arg=""
  [[ -n "${USE_LONG_ID:-}" ]] && long_arg="--long"
  local mm_all_hits="${out_mm}/mm_all_hits_${TASK_ID}.txt"
  : > "$mm_all_hits"

  run_mm_iteration() {
    local mm_out="${out_mm}/mm_iteration_${iteration}_${TASK_ID}.fasta"
    local mm_hit_ids="${out_mm}/mm_hits_${TASK_ID}.txt"
    local mm_fixed_out="${out_mm}/fixed_${TASK_ID}.fasta"
    local n_db n_fixed e_fixed e_db

    n_db=$(grep -c ">" "$mm_db_file" || true)
    n_fixed=$(grep -c ">" "$FIXED_FILE" || true)
    e_db=$(echo "scale=10; $E_VALUE * $n_db / $Z_SIZE" | bc)
    e_fixed=$(echo "scale=10; $E_VALUE * $n_fixed / $Z_SIZE" | bc)

    mkdir -p "${out_mm}/mm_db_${TASK_ID}" "${out_mm}/mm_fixed_${TASK_ID}" "${out_mm}/tmp_${TASK_ID}"

    $MMSEQS createdb "$FIXED_FILE" "${out_mm}/mm_fixed_${TASK_ID}/db" >&2
    $MMSEQS createdb "$mm_db_file" "${out_mm}/mm_db_${TASK_ID}/db" >&2

    if $MMSEQS search "${out_mm}/mm_fixed_${TASK_ID}/db" "${out_mm}/mm_db_${TASK_ID}/db" "${out_mm}/tmp_${TASK_ID}/db" "${out_mm}/tmp_fi_db_${TASK_ID}" --alignment-mode 2 --cov-mode 0 -c 0 -e "$e_fixed" --threads "$MMSEQS_CORES" --max-seqs "$MMSEQS_MAX_SEQS" >&2; then
      $MMSEQS convertalis "${out_mm}/mm_fixed_${TASK_ID}/db" "${out_mm}/mm_db_${TASK_ID}/db" "${out_mm}/tmp_${TASK_ID}/db" "${out_mm}/mm_hits_fi_db_${TASK_ID}.tsv" --format-mode 2 >&2
    else
      echo "Error: mmseqs search failed (Fixed -> DB)" >&2
      exit 1
    fi
    rm -rf "${out_mm}/tmp_${TASK_ID}"/*

    if $MMSEQS search "${out_mm}/mm_db_${TASK_ID}/db" "${out_mm}/mm_fixed_${TASK_ID}/db" "${out_mm}/tmp_${TASK_ID}/db" "${out_mm}/tmp_db_fi_${TASK_ID}" --alignment-mode 2 --cov-mode 0 -c 0 -e "$e_db" --threads "$MMSEQS_CORES" --max-seqs "$MMSEQS_MAX_SEQS" >&2; then
      $MMSEQS convertalis "${out_mm}/mm_db_${TASK_ID}/db" "${out_mm}/mm_fixed_${TASK_ID}/db" "${out_mm}/tmp_${TASK_ID}/db" "${out_mm}/mm_hits_db_fi_${TASK_ID}.tsv" --format-mode 2 >&2
    else
      echo "Error: mmseqs search failed (DB -> Fixed)" >&2
      exit 1
    fi
    rm -rf "${out_mm}/tmp_${TASK_ID}"*
    rm -rf "${out_mm}/mm_fixed_${TASK_ID}" "${out_mm}/mm_db_${TASK_ID}"

    if [[ "$REMOVE_TARGET" == "fixed" ]]; then
      # Fixed-side IDs: query col of fi_db (fixed=query), target col of db_fi (fixed=target)
      (awk '{print $1}' "${out_mm}/mm_hits_fi_db_${TASK_ID}.tsv"
       awk '{print $2}' "${out_mm}/mm_hits_db_fi_${TASK_ID}.tsv") | sort -u > "$mm_hit_ids"
    else
      # DB-side IDs: target col of fi_db, query col of db_fi
      (awk '{print $1}' "${out_mm}/mm_hits_db_fi_${TASK_ID}.tsv"
       awk '{print $2}' "${out_mm}/mm_hits_fi_db_${TASK_ID}.tsv") | sort -u > "$mm_hit_ids"
    fi
    cat "$mm_hit_ids" >> "$mm_all_hits"

    local line_count
    line_count=$(wc -l < "$mm_hit_ids")
    if [[ "$line_count" -gt 0 ]]; then
      if [[ "$REMOVE_TARGET" == "fixed" ]]; then
        remove_seqs_bash "$FIXED_FILE" "$mm_hit_ids" "$mm_fixed_out" $long_arg
        FIXED_FILE="$mm_fixed_out"
        check_fasta_nonempty "$mm_fixed_out" "MMseqs iteration $iteration"
      else
        remove_seqs_bash "$mm_db_file" "$mm_hit_ids" "$mm_out" $long_arg
        mm_db_file="$mm_out"
        check_fasta_nonempty "$mm_out" "MMseqs iteration $iteration"
      fi
    else
      echo "MMseqs iteration $iteration: no hits, skipping remove."
    fi
    MMSEQS_ITERATION_LINE_COUNT=$line_count
  }

  while true; do
    local line_count
    run_mm_iteration
    line_count=$MMSEQS_ITERATION_LINE_COUNT
    echo "MMseqs hits in iteration $iteration: $line_count"
    if [[ "$line_count" -eq 0 ]]; then
      echo "MMseqs stopping at iteration $iteration."
      break
    fi
    iteration=$((iteration + 1))
  done

  if [[ "$REMOVE_TARGET" == "db" ]]; then
    DB_FILE="$mm_db_file"
    check_fasta_nonempty "$DB_FILE" "MMseqs"
  else
    check_fasta_nonempty "$FIXED_FILE" "MMseqs"
  fi

  if [[ -z "$KEEP_INTERMEDIATES" ]]; then
    rm -f "${out_mm}"/mm_hits_fi_db_${TASK_ID}.tsv "${out_mm}"/mm_hits_db_fi_${TASK_ID}.tsv "${out_mm}"/mm_hits_${TASK_ID}.txt
  fi
  echo "MMSEQS took $((SECONDS - start)) seconds"
}

# --- BLAST ---
run_blast() {
  local start=$SECONDS
  local out_blast="${OUT_DIR}/blast_${OUT_SUFFIX}"
  mkdir -p "$out_blast"
  local blast_out="${out_blast}/blast_${TASK_ID}.fasta"
  local blast_hit_ids="${out_blast}/blast_hits_${TASK_ID}.txt"
  mkdir -p "${out_blast}/tmp_${TASK_ID}" "${out_blast}/b_fixed_${TASK_ID}" "${out_blast}/b_db_${TASK_ID}"

  $BLAST_DIR/makeblastdb -in "$FIXED_FILE" -dbtype prot -out "${out_blast}/b_fixed_${TASK_ID}/db" >&2
  $BLAST_DIR/makeblastdb -in "$DB_FILE" -dbtype prot -out "${out_blast}/b_db_${TASK_ID}/db" >&2

  $BLAST_DIR/blastp -query "$FIXED_FILE" -db "${out_blast}/b_db_${TASK_ID}/db" \
    -out "${out_blast}/blast_hits_fi_db_${TASK_ID}.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" -evalue "$E_VALUE" -dbsize "$BLAST_DBSIZE" -max_target_seqs "$BLAST_MAX_TARGET_SEQS" -num_threads "$BLAST_CORES" >&2
  $BLAST_DIR/blastp -query "$DB_FILE" -db "${out_blast}/b_fixed_${TASK_ID}/db" \
    -out "${out_blast}/blast_hits_db_fi_${TASK_ID}.tsv" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" -evalue "$E_VALUE" -dbsize "$BLAST_DBSIZE" -max_target_seqs "$BLAST_MAX_TARGET_SEQS" -num_threads "$BLAST_CORES" >&2

  rm -rf "${out_blast}/tmp_${TASK_ID}"* "${out_blast}/b_fixed_${TASK_ID}" "${out_blast}/b_db_${TASK_ID}"

  (awk '{print $1}' "${out_blast}/blast_hits_db_fi_${TASK_ID}.tsv"
   awk '{print $2}' "${out_blast}/blast_hits_fi_db_${TASK_ID}.tsv") | sort -u > "$blast_hit_ids"

  if [[ "$REMOVE_TARGET" == "fixed" ]]; then
    remove_seqs_bash "$FIXED_FILE" "$blast_hit_ids" "${out_blast}/fixed_blast_${TASK_ID}.fasta"
    FIXED_FILE="${out_blast}/fixed_blast_${TASK_ID}.fasta"
    check_fasta_nonempty "$FIXED_FILE" "BLAST"
  else
    remove_seqs_bash "$DB_FILE" "$blast_hit_ids" "$blast_out"
    DB_FILE="$blast_out"
    check_fasta_nonempty "$blast_out" "BLAST"
  fi

  if [[ -z "$KEEP_INTERMEDIATES" ]]; then
    rm -f "${out_blast}/blast_hits_fi_db_${TASK_ID}.tsv" "${out_blast}/blast_hits_db_fi_${TASK_ID}.tsv"
  fi
  echo "BLAST took $((SECONDS - start)) seconds"
}

# --- SW (Smith-Waterman) ---
run_sw() {
  local start=$SECONDS
  local out_sw="${OUT_DIR}/sw_${OUT_SUFFIX}"
  mkdir -p "$out_sw"
  local sw_out="${out_sw}/sw_${TASK_ID}.fasta"
  local sw_hit_ids="${out_sw}/sw_hits_${TASK_ID}.txt"

  $FASTA_DIR/ssearch36 -m 8 -T "$SW_CORES" -E "$E_VALUE" -Z "$Z_SIZE" \
    "$FIXED_FILE" "$DB_FILE" > "${out_sw}/sw_hits_fi_db_${TASK_ID}.tsv"
  $FASTA_DIR/ssearch36 -m 8 -T "$SW_CORES" -E "$E_VALUE" -Z "$Z_SIZE" \
    "$DB_FILE" "$FIXED_FILE" > "${out_sw}/sw_hits_db_fi_${TASK_ID}.tsv"

  (awk '{print $1}' "${out_sw}/sw_hits_db_fi_${TASK_ID}.tsv"
   awk '{print $2}' "${out_sw}/sw_hits_fi_db_${TASK_ID}.tsv") | sort -u > "$sw_hit_ids"

  if [[ "$REMOVE_TARGET" == "fixed" ]]; then
    remove_seqs_bash "$FIXED_FILE" "$sw_hit_ids" "${out_sw}/fixed_sw_${TASK_ID}.fasta"
    FIXED_FILE="${out_sw}/fixed_sw_${TASK_ID}.fasta"
    check_fasta_nonempty "$FIXED_FILE" "Smith-Waterman"
  else
    remove_seqs_bash "$DB_FILE" "$sw_hit_ids" "$sw_out"
    DB_FILE="$sw_out"
    check_fasta_nonempty "$sw_out" "Smith-Waterman"
  fi

  if [[ -z "$KEEP_INTERMEDIATES" ]]; then
    rm -f "${out_sw}/sw_hits_fi_db_${TASK_ID}.tsv" "${out_sw}/sw_hits_db_fi_${TASK_ID}.tsv"
  fi
  echo "SW took $((SECONDS - start)) seconds"
}

# --- Dispatch by --order (characters validated before config load) ---
for (( i=0; i<${#ORDER}; i++ )); do
  case "${ORDER:i:1}" in
    p|P) run_phmmer ;;
    m|M) run_mmseqs ;;
    b|B) run_blast ;;
    s|S) run_sw ;;
    *) echo "Unknown tool in --order: ${ORDER:i:1}" >&2; exit 1 ;;
  esac
done

echo "Completed. Removed from $REMOVE_TARGET  DB_FILE=$DB_FILE  FIXED_FILE=$FIXED_FILE"
