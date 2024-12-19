#!/bin/bash

mkdir -p build 

cmake -S . \
      -B build \
      -G Ninja \
      -DCMAKE_CXX_COMPILER=clang++-18 \
      -DCMAKE_C_COMPILER=clang-18 \
      -DTYPES="FLOAT,FAST_FIXED(13,7),FAST_FIXED(32,5),DOUBLE" \
      -DSIZES="S(36,84),S(10,10),S(42,1337)" \
      -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping -O2" \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=1
#-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping" \
cmake --build build --parallel $(nproc)
cp -f build/compile_commands.json .

./build/fluid --p-type=FLOAT --v-type=FAST_FIXED\(13,7\) --v-flow-type=FAST_FIXED\(32,5\) --field="fields/example1"
