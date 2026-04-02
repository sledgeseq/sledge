# Sledge

Sledge splits protein sequence databases into train / test / (optional) validation sets (`sledge_splitter`), scores and filters sequence pairs with a pHMMER-based engine (`phmmer_filter`), and orchestrates multi-tool filtering pipelines (`sledge_filter`: pHMMER, MMseqs2, BLAST, FASTA `ssearch36`). Sledge is released under the MIT License (see [License](#license)).

## Citation

To cite Sledge, please use:

> [SLEDGE_CITATION_PLACEHOLDER]
>
> [SLEDGE_PAPER_LINK_PLACEHOLDER]

---

## Acknowledgements

We would like to thank the developers of [HMMER3](http://hmmer.org) <sup><a href="#ref-eddy-2011">1</a></sup> for the EASEL tools and base pHMMER pipeline that have been customized for `sledge_splitter` and `sledge_filter`.

---

## OS support matrix

| OS | Build support | `sledge_filter` support | External tools support | Notes |
|----|---------------|-------------------------|------------------------|-------|
| Linux (x86_64) | Official | Official | Official | Primary supported platform. |
| macOS (Apple Silicon / Intel) | Official | Official | Official | Use Homebrew/system packages for dependencies. |
| Windows (WSL2) | Official (via Linux in WSL2) | Official (via Bash in WSL2) | Official (via Linux in WSL2) | Recommended Windows path. |
| Windows (native) | Planned | Planned | Planned | Not currently first-class; requires separate installers and shell adaptations. |

---

## Requirements

**Core build tools (all supported platforms)**  
- C compiler (`gcc`/`clang`)
- GNU `make`
- `ar`, `ranlib`
- Bash (for `sledge_filter` and install helpers)
- `pthread` and math library (`-lm`)

**Runtime dependencies for `sledge_filter`**  
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (`makeblastdb`, `blastp`)
- [FASTA package](https://github.com/wrpearson/fasta36) (`ssearch36`)

Refer to the [external tools](#external-tools) section for installation through `make`.

---

## Installation

### 1) Clone and build core binaries

```bash
git clone https://github.com/YOUR_ORG/sledge.git
cd sledge/ # adjust if your repo root differs
./configure   # detects CPU (SSE on x86_64, NEON on aarch64/arm64); writes build-config.mk
make
```

**Shell (bash vs zsh):** It does not matter whether you use bash or zsh as your login shell. `./configure` has a `#!/usr/bin/env bash` shebang, so it always runs under bash. The same commands work in Terminal with zsh.

**If you see “Permission denied” when running `./configure`:** the file may not be executable in your clone (e.g. checkout without the executable bit). Fix with `chmod +x configure`, or run `bash configure` (no execute bit needed).

Re-run `./configure` after changing machines or if you want to force a backend (`./configure --with-impl=sse` or `--with-impl=neon`). `make distclean` removes `build-config.mk` as well as build products.

Executables land in `bin/`:

| Binary | Role |
|--------|------|
| `sledge_splitter` | Profmark-style split of one FASTA DB into train / test / val |
| `phmmer_filter` | Pairwise sequence similarity filter using HMMER |
| `sledge_filter` | Shell pipeline: configurable order of p / m / b / s steps |

Add `bin` to your `PATH` or call tools with full paths:

```bash
export PATH="/path/to/sledge/bin:$PATH"
```

### 2) Install external tools (MMseqs2, BLAST+, FASTA/ssearch36)

Install external tools through the Makefile:

```bash
make install-external
```

### External tools

This installs MMseqs2, BLAST+, and FASTA tools into `sledge/external_tools` by default.

**Linux** (`install/linux/install_external_linux.sh`) selects downloads from `uname -m`: **x86_64** uses MMseqs2 `sse2` binaries and NCBI **x64-linux** BLAST; **aarch64**/**arm64** uses MMseqs2 **arm64** and NCBI **aarch64-linux** BLAST, and builds FASTA36 with non-SSE makefiles. Override MMseqs flavor with `MMSEQS_ARCH` (e.g. `avx2` on capable x86_64) or NCBI tarball name with `BLAST_PLATFORM` if needed.

If external tools are already installed, point `sledge_filter` to them in your config (see `install/test_filter.config`):

| Config key | Set this to |
|------------|-------------|
| `MMSEQS` | Full path to `mmseqs` executable |
| `BLAST_DIR` | Directory containing `makeblastdb` and `blastp` |
| `FASTA_DIR` | Directory containing `ssearch36` |

`sledge_filter` can also run on a subset of tools. Only the tools used in `--order` need to be installed.

`make install-external` uses the cross-platform installer dispatcher under `install/`.
You can still call the dispatcher directly if needed:

```bash
./install/install_external.sh /path/to/sledge
```

The [fasta36](https://github.com/wrpearson/fasta36) checkout used by the installer needs compatibility patches to build with current GCC and Clang; see [FASTA36 GCC patches](#fasta36-gcc-patches) near the end of this README.

**macOS:** `install/macos/install_external_macos.sh` installs MMseqs2 and BLAST+ via **Homebrew** (`brew`, native for Intel vs Apple Silicon) when available, and otherwise falls back to direct upstream tarball downloads (`mmseqs-osx-universal`, `ncbi-blast-<version>+-x64-macosx` by default). FASTA36 is built from source, choosing FASTA makefiles by **`uname -m`** (ARM Macs no longer prefer x86-only makefiles first). Optional overrides: `MMSEQS_TAG`, `MMSEQS_MACOS_ASSET`, `BLAST_VERSION`, `BLAST_PLATFORM`.

On Windows, run the PowerShell entrypoint (WSL2-backed flow):

```powershell
powershell -ExecutionPolicy Bypass -File .\install\windows\install_external_windows.ps1 -SledgeDir C:\path\to\sledge
```

### 3) Run installation tests

After building binaries and configuring/installing external tools:

```bash
make test-install
```

`make test-install` auto-detects OS:
- Linux/macOS: runs `install/test_installation.sh`
- Windows: routes to `install/windows/test_installation_windows.ps1` (WSL2-backed)

Direct Windows entrypoint is still available if needed:

```powershell
powershell -ExecutionPolicy Bypass -File .\install\windows\test_installation_windows.ps1 -SledgeDir C:\path\to\sledge
```

---

## Quick start

**Splitter — one input DB, written train/test/val under the current directory**

```bash
sledge_splitter --dbblock 100 --test_limit 20 --val_limit 10 -o stats --output_dir results seqDB.fasta
```

**Multi-tool pipeline — order `pmbs` (pHMMER → MMseqs2 → BLAST → Smith–Waterman)**

```bash
sledge_filter --order pmbs --config pipeline.config --fixed-file test.fasta --db-file train_candidates.fasta
```

---

## `sledge_splitter`

### Required arguments

| Argument | Description |
|----------|-------------|
| `<seqdb>` | One input protein sequence file (positional; must be last). |

### Essential options (typical runs)

| Option | Default | Description |
|--------|---------|---------------|
| `-o <prefix>` | `-` | Prefix for stats / summary output files (e.g. `stats` → `stats_0.txt`). |
| `-Z <n>` | *inferred from database* | Effective database size for E-value calculation. |
| `--cpu <n>` | `1` | Worker threads. |
| `--dbblock <n>` | `1000` | Number of sequences in database. |
| `--test_limit <n>` | `75000` | Stop once at least this many sequences are assigned to **test**. |
| `--val_limit <n>` | `10000` | Stop once at least this many are assigned to **validation** (unless disabled). |
| `--init_chunk <n>` | `10` | Sequences considered for assignment at one go |
| `--seed <n>` | `42` | RNG seed (`0` = one-time random seed). |
| `--suppress` | off | Disable progress bar. |
| `--task_id <id>` | `0` | Suffix for output files (`*_0.fasta`, etc.). |

### Output paths

| Option | Description |
|--------|-------------|
| `--output_dir <dir>` | Write `train`, `test`, `val`, `discard` under `<dir>` (with default basenames unless overridden). |
| `--train_path`, `--test_path`, `--val_path`, `--discard_path` | Basename prefixes for the four FASTA streams (defaults: `train`, `test`, `val`, `discard`). |

### Similarity / scoring (pHMMER-style)

| Option | Description |
|--------|-------------|
| `-E <x>` | Treat hits with E-value ≤ `x` as significant. |
| `--plow`, `--phigh` | PID window for accepting a sequence into the train/test/val set (see `-h` for interaction with the algorithm). |

### Resume / debug

| Option | Description |
|--------|-------------|
| `--load_tr`, `--load_te`, `--load_val` | Resume from existing train/test/val FASTA (`-` = none). |
| `--halt <n>` | Process at most `n` sequences (debug). |
| `--resume <n>` | Start after sequence index `n`. |
| `--freq <n>` | Progress update frequency. |

For complete set of options: run `sledge_splitter -h`.

---

## `sledge_filter` (pipeline script)

### Required CLI arguments

| Option | Description |
|--------|-------------|
| `--order <string>` | Tool order: `p` = pHMMER (`phmmer_filter`), `m` = MMseqs2, `b` = BLAST+, `s` = Smith–Waterman (`ssearch36`). Example: `pmbs`, `spm`. |
| `--config <file>` | Bash-sourced config (see below). |
| `--fixed-file <fasta>` | “Fixed” side (e.g. reference or test set). |
| `--db-file <fasta>` | Database to filter against the fixed set. |

### Optional CLI

| Option | Description |
|--------|-------------|
| `--out-suffix <name>` | Suffix for per-tool output dirs (`phmmer_<suffix>`, …); defaults to `TASK_ID` from config. |

### Config file (required variables)

`OUT_DIR` must be set (base directory for outputs). Typically also set:

| Variable | Role |
|----------|------|
| `PHMMER_FILTER` | Path to `phmmer_filter` binary. |
| `MMSEQS` | MMseqs2 executable. |
| `BLAST_DIR` | Directory containing `makeblastdb` and `blastp`. |
| `FASTA_DIR` | Directory containing `ssearch36`. |
| `TASK_ID` | Integer or `SLURM` (uses `SLURM_ARRAY_TASK_ID`). |
| `REMOVE_TARGET` | `db` or `fixed` — which side sequences are removed from after hits. |
| `E_VALUE` | E-value threshold (shared across tools where applicable). |
| `Z_SIZE` | Positive integer; used for pHMMER / SW calibration. |
| `SLEDGE_CORES`, `MMSEQS_CORES`, `BLAST_CORES`, `SW_CORES` | Thread counts for each tool. |
| `BLAST_DBSIZE`, `BLAST_MAX_TARGET_SEQS`, `MMSEQS_MAX_SEQS` | BLAST / MMseqs2 tuning. |

Optional: `KEEP_INTERMEDIATES`, `USE_LONG_ID`, etc. (see script header in `src/sledge_filter.sh`).

**Minimal example config**

```bash
OUT_DIR="./filter_run_out"
PHMMER_FILTER="/path/to/sledge_minimal/bin/phmmer_filter"
MMSEQS="mmseqs"
BLAST_DIR="/path/to/ncbi-blast/bin"
FASTA_DIR="/path/to/fasta/bin"
TASK_ID=0
```

---

## `phmmer_filter`

If you are not interested in using a suite of tools to filter and want to use just the optimized PHMMER filter for different functionalities, refer to the documentation below.

### Required arguments

| Argument | Description |
|----------|-------------|
| `<qdb>` | Query sequence database (first positional). |
| `<tdb>` | Target sequence database (second positional). |

One of `qdb` or `tdb` may be `-` (stdin), not both.

### Essential options

| Option | Default | Description |
|--------|---------|-------------|
| `-o <prefix>` | `-` | Output prefix for result files. |
| `-Z <n>` | *inferred from database by default* | Database size for E-value calibration. |
| `-E <x>` | `10.0` | Reporting E-value threshold. |
| `--cpu <n>` | `1` | Threads. |
| `--qsize <n>` | `1` | Queries per thread per batch. |
| `--qblock <n>` | `10` | Query block size. |
| `--tblock <n>` | `1000` | Target block size. |
| `--format <n>` | `1` | Output format (see below). |
| `--all_hits` | off | Do not stop early on first failing comparison (when applicable). |
| `--task_id <id>` | `0` | Shard id in output filenames. |
| `--seed <n>` | `42` | RNG seed. |
| `--suppress` | off | Disable progress bar. |

### Output formats (`--format`)

| Value | Meaning |
|------|---------|
| `0` | Per-sequence ACCEPT/REJECT string. |
| `1` | Full information for rejected hits (default). |
| `2` | IDs of accepted sequences. |
| `3` | IDs of rejected queries. |
| `4` | IDs of rejected targets (e.g. with `--all_hits`). |

### Input formats

`--qformat` / `--tformat` force Easel format names (e.g. `fasta`) and skip autodetection.

---

## Example use cases

1. **Splitter — train/val/test from a single database**  
   Build disjoint splits for benchmarking: tune `--test_limit`, `--val_limit`, `-Z`, and `--cpu` for large DBs.

   ```bash
   sledge_splitter -Z 1e6 --cpu 16 --test_limit 5000 --val_limit 1000 -o run stats_db.fasta
   ```

2. **Filter — large dissimilar training set given a fixed test set**  
   Run `sledge_filter` with `REMOVE_TARGET=db` and strict `E_VALUE` / `Z_SIZE` so the “db” FASTA loses sequences similar to the test set (after pHMMER / MMseqs / BLAST / SW as ordered).

3. **Filter — remove overlap between any two FASTA sets**  
   Point `--fixed-file` and `--db-file` at the two pools; choose `REMOVE_TARGET` to drop hits from either side.

4. **Splitter / filter — out-of-distribution (OOD) evaluation**  
   Use the splitter to hold out a test set; use `sledge_filter` or `phmmer_filter` to strip training candidates that match the test set, giving a cleaner OOD evaluation. Alternatively, prune the test set with respect to the train set to get a subset of OOD sequences. You can also deduplicate a set of sequences by removing ones that have high sequence similarity to each other within a set (this would involve filtering a set with itself).

5. **Strictness, order, and tool subset**  
   - **Strictness:** lower `E_VALUE`, appropriate `Z_SIZE`, and `BLAST_DBSIZE` / MMseqs e-value scaling in the script.  
   - **Order:** e.g. `p` only for fast homology removal; `pmbs` for a long cascade.  
   - **Speed** MMSeqs2 is the fastest of the tools, followed by BLAST and PHMMER (our modified implementation) in similar order of magnitude. Smith-Waterman based filtering by ssearch is the slowest part of this process. User can choose the tools and order as required.
   - **Subset:** omit letters (e.g. `pm` skips BLAST and SW).

6. **Stratification / reporting**  
   Use `phmmer_filter` `--format` modes and `-E` / `-Z` to export accept/reject IDs or full hit tables; combine with downstream scripts to bin by PID or E-value from the tabular output. Our customized pHMMER filter allows reporting results in different forms. We have not modified the source code of MMSeqs2, BLAST and ssearch36.

---

## Further help

```bash
sledge_splitter -h
phmmer_filter -h
```

For sledge_filter pipeline behavior and defaults, read the comments at the top of `src/sledge_filter.sh`.

---

## FASTA36 GCC patches

Upstream [wrpearson/fasta36](https://github.com/wrpearson/fasta36) predates strict ISO C defaults and C23 keywords in recent compilers. Sledge ships a unified diff so the cloned sources compile cleanly on typical Linux and macOS toolchains.

### Patch file and tooling

| Item | Location / note |
|------|-----------------|
| Patch | `install/patches/fasta36-gcc-prototypes.patch` |
| Applied by | Linux/macOS installers after `git clone` (see `install/linux/install_external_linux.sh`, `install/macos/install_external_macos.sh`) |
| System dependency | The **`patch`** command must be on `PATH` (install the `patch` package on minimal images if needed). The patch *file* is part of this repo; only the `patch` binary is external. |

### Scope (13 files under `src/`)

Patches are **cumulative**: partial fixes (e.g. only mmap prototypes) are not enough. Grouped by theme:

- **Headers / mmap:** `altlib.h`, `mm_file.h`, `mmgetaa.c` — correct prototypes and function-pointer types for library and mmap readers.
- **C23 / syntax:** `pssm_asn_subs.c` — avoid the `bool` keyword clash; full `parse_pssm_*` prototypes. `lsim4.h` — `#include <stdbool.h>` instead of `typedef int bool`.
- **Tool front-ends:** `map_db.c` — typed function pointers for `get_entry` / `get_ent_arr`. `initfa.c` — ANSI `sortbest`, safer `get_lambda` allocation. `compacc2e.c` — bounded `calloc` for annotation tables; on UNIX, `mktemp` replaced with `mkstemp` for the temp library DB path.
- **FastA/FastX drops:** `dropfx2.c`, `dropfz3.c` — K&R-style forward declarations, `lx_band` / `ckalloc`, `kssort`, `global` / `small_global` / `global_up` / `global_down` (including `const` sequence pointers where needed), `fatal`. `dropfs2.c`, `dropff2.c`, `dropnfa.c` — ANSI prototypes for shell-sort helpers (`kssort`, `kpsort`, `krsort`) and `savemax` where applicable.

### Branch and maintenance

- Installers default to **`FASTA36_REF=master`** (see `install/install_external.sh`). The patch is generated against **`master`**. If GitHub’s default branch or your chosen tag diverges, line numbers may not match—regenerate the diff from a clean checkout of that revision, or use `master`.
- If upstream edits the same lines on `master`, refresh the patch from a fresh clone and re-run the install build.

### Hosts without `patch`

Nothing in the repo vendors a `patch` binary today. For air-gapped or minimal hosts, you could ship a static GNU `patch` (e.g. under `install/bin/`) and invoke it from the installers.

---

## References

1. <a id="ref-eddy-2011"></a>Eddy SR (2011) Accelerated Profile HMM Searches. *PLOS Computational Biology* 7(10): e1002195. https://doi.org/10.1371/journal.pcbi.1002195

---

## License

This project is licensed under the MIT License. See `LICENSE` for the full text.
