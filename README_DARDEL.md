# Building and Running HiCapTools on Dardel

This guide builds the optimized HiCapTools version with the restriction-site
boundary fix. The build uses Dardel's Boost headers and the bundled BamTools
library.

## 1. Get the source

Clone the repository and enter it:

```bash
git clone --depth 1 https://github.com/sahlenlab/HiCapTools.git
cd HiCapTools
```

The repository contains `CMakeLists.txt`, `src`, `include`, `scripts`,
`bamtools`, and `bin`. Boost is intentionally not bundled; the build uses
Dardel's Boost module.

## 2. Load the build environment

Start from a fresh login shell when possible, then load:

```bash
ml PDC
ml PrgEnv-gnu
ml boost
ml cmake/3.31.3
ml
```

Check the tools and Boost installation:

```bash
cmake --version
CC --version
echo "$EBROOTBOOST"
```

`EBROOTBOOST` must not be empty.

## 3. Compile in Release mode

The supplied script configures a release build, uses the bundled Linux
BamTools library, and finds Boost from the loaded Dardel module:

```bash
chmod +x scripts/build_dardel.sh
./scripts/build_dardel.sh
```

The executable is written to:

```text
bin/HiCapTools
```

The final build output should contain:

```text
[100%] Built target HiCapTools
```

## 4. Verify linked libraries

```bash
file bin/HiCapTools
ldd bin/HiCapTools
```

The executable must be an ELF 64-bit Linux executable, and `ldd` must not show
any `not found` entries. BamTools should resolve inside this installation:

```bash
ldd bin/HiCapTools | grep bamtools
```

## 5. Prepare the configuration

Absolute paths are recommended for all cluster inputs, especially:

- Experiment file
- BAM files and their `.bai` indexes
- Probe and negative-control probe files
- Digested genome file
- Transcript and SNV files
- Negative-control regions
- ENCODE blacklist

When an existing digested genome file is configured, `hg38.fa` is not needed.

To enable blacklist filtering and integrated output, include:

```text
Blacklist File=/absolute/path/ENCODE.BlackListedRegions.v2.hg38.txt
Generate Integrated Interactions=Yes
```

Confirm that the run log contains `Blacklist regions loaded:`. If the file is
not accessible, HiCapTools warns and continues without blacklist filtering.

The main and negative-control integrated files contain a
`MergedInteractorID` column. Distal interacting intervals on the same
chromosome are assigned to the same merged region when they overlap or are at
most 150 bases apart; this is applied transitively. Every original interaction
remains on its own row. Probe-probe interactions retain their original
`InteractorID` in this column.

## 6. Submit a job

The supplied job script defaults to `HiCapTools` in the same `bin` directory
as the script. To use another installation, set its path before submitting:

```bash
export HICAPTOOLS_EXE=/absolute/path/to/HiCapTools/bin/HiCapTools
```

Move to the directory where output should be written, then submit chromosome
22:

```bash
cd /path/to/interactionCalls
sbatch -A YOUR_DARDEL_PROJECT \
    /absolute/path/to/HiCapTools/bin/runHiCapTools_dardel.sbatch \
    chr22 /absolute/path/to/config/configFile.txt
```

Use `All` instead of `chr22` to process all chromosomes. If the submission
directory contains `config/configFile.txt`, the config argument may be omitted:

```bash
sbatch -A YOUR_DARDEL_PROJECT \
    /absolute/path/to/HiCapTools/bin/runHiCapTools_dardel.sbatch All
```

Replace `YOUR_DARDEL_PROJECT` with the allocation that should be charged.

Monitor the job with:

```bash
squeue -u "$USER"
tail -f HiCapTools-JOBID.out
sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

A successful log ends with `Execution Complete` and `Finished:`.

## Troubleshooting

### Boost is not detected

Reload the environment and confirm `EBROOTBOOST`:

```bash
ml PDC
ml PrgEnv-gnu
ml boost
echo "$EBROOTBOOST"
```

### CMake policy error

Use the documented CMake module:

```bash
ml cmake/3.31.3
```

The current source requires CMake 3.10 or newer.

### A linked library is missing at runtime

Load the same PDC/GNU environment used for the build and inspect:

```bash
ldd bin/HiCapTools | grep "not found"
```
