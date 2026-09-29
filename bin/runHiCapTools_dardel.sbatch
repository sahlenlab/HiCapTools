#!/bin/bash -l
#SBATCH -J HiCapTools
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -o HiCapTools-%j.out
#SBATCH -e HiCapTools-%j.err

set -euo pipefail

chromosome=${1:-chr22}
config_file=${2:-${SLURM_SUBMIT_DIR}/config/configFile.txt}
mode=${HICAP_MODE:-PrintProximities}
proximity_type=${HICAP_PROXIMITY_TYPE:-Both}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
executable=${HICAPTOOLS_EXE:-${script_dir}/HiCapTools}

ml PDC
ml PrgEnv-gnu
ml boost

cd "${SLURM_SUBMIT_DIR}"

if [[ ! -x "${executable}" ]]; then
    echo "HiCapTools executable not found or not executable: ${executable}" >&2
    exit 1
fi

if [[ ! -f "${config_file}" ]]; then
    echo "Config file not found: ${config_file}" >&2
    exit 1
fi

# RPATH is set during the Dardel build; this is a fallback for older loaders.
project_root=$(cd "$(dirname "${executable}")/.." && pwd)
export LD_LIBRARY_PATH="${project_root}/bamtools${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=1

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Executable: ${executable}"
echo "Mode: ${mode}"
echo "Chromosome: ${chromosome}"
echo "Config: ${config_file}"

if [[ "${mode}" == "ComputeStatsOnly" ]]; then
    srun --ntasks=1 "${executable}" \
        ProximityDetector -m ComputeStatsOnly -config "${config_file}"
else
    srun --ntasks=1 "${executable}" \
        ProximityDetector -m PrintProximities -p "${proximity_type}" \
        -c "${chromosome}" -config "${config_file}"
fi

echo "Finished: $(date)"
