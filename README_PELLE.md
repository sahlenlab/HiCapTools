# Building and Running HiCapTools on Pelle

This guide builds the optimized HiCapTools version with the restriction-site
boundary fix. It uses Pelle's GCC 13.3 toolchain and Boost headers together
with the bundled BamTools library.

## 1. Get the source

Clone the repository and enter it:

```bash
git clone https://github.com/sahlenlab/HiCapTools.git
cd HiCapTools
```

The repository contains `CMakeLists.txt`, `src`, `include`, `bamtools`, and
`bin`. Boost is intentionally not bundled; the build uses Pelle's Boost
module.

## 2. Load the build environment

If Conda is active, deactivate it before loading compiler modules:

```bash
conda deactivate
```

If this reports that Conda is not active, continue. Load one consistent module
family:

```bash
module purge
module load GCC/13.3.0
module load CMake/3.31.8-GCCcore-13.3.0
module load Boost/1.85.0-GCC-13.3.0
module list
```

Verify the environment:

```bash
g++ --version
cmake --version
echo "$EBROOTBOOST"
```

The expected compiler and CMake versions are GCC 13.3.0 and CMake 3.31.8.
`EBROOTBOOST` must not be empty.

## 3. Compile in Release mode

Configure a Pelle-specific build directory:

```bash
cmake -S . -B build-pelle \
    -DCMAKE_BUILD_TYPE=Release \
    -DHICAPTOOLS_USE_SYSTEM_BAMTOOLS=OFF \
    -DBOOST_INCLUDE_DIR="$EBROOTBOOST/include"
```

Compile:

```bash
cmake --build build-pelle --parallel 4
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
ldd bin/HiCapTools | grep "not found"
```

The executable must be an ELF 64-bit Linux executable. The final command must
produce no output. Confirm that BamTools resolves inside this installation:

```bash
ldd bin/HiCapTools | grep bamtools
```

## 5. Prepare the configuration

Do not reuse Dardel paths beginning with `/cfs/klemming/`. Update the config and
experiment file to use paths that exist on Pelle. Absolute paths are
recommended for:

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
Blacklist File=/absolute/pelle/path/ENCODE.BlackListedRegions.v2.hg38.txt
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

Set the actual executable and job-script paths while inside the installation:

```bash
export HICAPTOOLS_EXE="$(realpath bin/HiCapTools)"
export HICAPTOOLS_SBATCH="$(realpath bin/runHiCapTools_pelle.sbatch)"
```

Move to the directory where output should be written, then submit chromosome
22:

```bash
cd /path/to/interactionCalls
sbatch -A YOUR_PELLE_PROJECT \
    "$HICAPTOOLS_SBATCH" chr22 /absolute/path/to/config/configFile.txt
```

Use `All` instead of `chr22` to process all chromosomes. If the submission
directory contains `config/configFile.txt`, the config argument may be omitted:

```bash
sbatch -A YOUR_PELLE_PROJECT "$HICAPTOOLS_SBATCH" All
```

The supplied Pelle script loads GCC 13.3 and its matching zlib module at
runtime. Replace `YOUR_PELLE_PROJECT` with an active Pelle compute account.
List valid accounts with:

```bash
sacctmgr show assoc user="$USER" format=Account,Partition
```

Monitor the job with:

```bash
squeue -u "$USER"
tail -f HiCapTools-JOBID.out
sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

A successful log ends with `Execution Complete` and `Finished:`.

## Sharing with project members

For a project-owned folder that all group members may modify:

```bash
chmod -R g+rwX,o-rwx /path/to/shared/folder
find /path/to/shared/folder -type d -exec chmod g+s {} +
find /path/to/shared/folder -type d \
    -exec setfacl -d -m u::rwx,g::rwx,o::--- {} +
```

Add the following to Slurm scripts that create shared outputs:

```bash
umask 0007
```

## Troubleshooting

### Boost is not detected

```bash
module load Boost/1.85.0-GCC-13.3.0
echo "$EBROOTBOOST"
```

Then rerun the CMake configure command.

### A linked library is missing

Reload the same modules used for compilation:

```bash
module purge
module load GCC/13.3.0
module load CMake/3.31.8-GCCcore-13.3.0
module load Boost/1.85.0-GCC-13.3.0
ldd bin/HiCapTools | grep "not found"
```

### The job cannot find HiCapTools

Set the executable explicitly before `sbatch`:

```bash
export HICAPTOOLS_EXE=/absolute/path/to/HiCapTools/bin/HiCapTools
```
