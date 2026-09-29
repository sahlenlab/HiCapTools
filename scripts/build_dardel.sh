#!/bin/bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "${script_dir}/.." && pwd)
build_dir="${project_root}/build-dardel"
build_jobs=${SLURM_CPUS_PER_TASK:-4}

if ! command -v cmake >/dev/null 2>&1; then
    echo "cmake is not available. Load the Dardel CMake module first." >&2
    exit 1
fi

boost_root=${BOOST_HOME:-${BOOST_ROOT:-${EBROOTBOOST:-}}}
if [[ -z "${boost_root}" || ! -f "${boost_root}/include/boost/accumulators/accumulators.hpp" ]]; then
    echo "Dardel Boost headers are not loaded. Run 'ml PDC', 'ml PrgEnv-gnu', and 'ml boost' first." >&2
    exit 1
fi

cmake -S "${project_root}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DHICAPTOOLS_USE_SYSTEM_BAMTOOLS=OFF \
    -DBOOST_INCLUDE_DIR="${boost_root}/include"
cmake --build "${build_dir}" --parallel "${build_jobs}"

echo
echo "Built executable: ${project_root}/bin/HiCapTools"
file "${project_root}/bin/HiCapTools"
ldd "${project_root}/bin/HiCapTools"
