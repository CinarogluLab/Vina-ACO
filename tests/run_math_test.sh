#!/usr/bin/env bash
# Build and run the ACO math invariants (tests/test_aco_math.cpp).
set -euo pipefail
cd "$(dirname "$0")/.."
./build.sh -DVINA_ACO_TESTS=ON
cmake --build build --target test_aco_math
./build/test_aco_math
