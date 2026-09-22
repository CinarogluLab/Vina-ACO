#!/bin/bash
# One-command build of Vina-ACO -> bin/vina_aco
# Uses the Boost that CMake finds. To pick a specific installation, set BOOST_ROOT
# (for example BOOST_ROOT=$CONDA_PREFIX inside the conda environment of environment.yml).
set -euo pipefail
cd "$(dirname "$0")"
BOOST_ROOT=${BOOST_ROOT:-${CONDA_PREFIX:-}}
extra=()
if [ -n "$BOOST_ROOT" ]; then extra=(-DBOOST_ROOT="$BOOST_ROOT" -DBoost_NO_SYSTEM_PATHS=ON); fi
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release "${extra[@]}" "$@"
cmake --build build -j"$(nproc)"
echo "Built: $(pwd)/bin/vina_aco"
