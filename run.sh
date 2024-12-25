#!/bin/bash

mkdir -p build 

cmake -S . \
      -B build \
      -G Ninja \
      -DCMAKE_CXX_COMPILER=clang++-18 \
      -DCMAKE_C_COMPILER=clang-18 \
      -DTYPES="DOUBLE" \
      -DSIZES="S(36,84)" \
      -DCMAKE_CXX_FLAGS="-Wall -Wextra -pedantic -fopenmp=libomp" \
      -DCMAKE_BUILD_TYPE=release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=1 

#-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping -O2" \
cmake --build build --parallel $(nproc) || exit 1
cp -f build/compile_commands.json .

./build/fluid --p-type="DOUBLE" --v-type="DOUBLE" --v-flow-type="DOUBLE" --field="fields/example1" --j=12
