#!/bin/bash

mkdir -p build 

cmake -S . \
      -B build \
      -G Ninja \
      -DCMAKE_CXX_COMPILER=clang++ \
      -DCMAKE_C_COMPILER=clang \
      -DTYPES="FLOAT,FAST_FIXED(13,7),FIXED(32,5),DOUBLE" \
      -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping" \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=1 
#-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping" \
cmake --build build --parallel $(nproc)
cp -f build/compile_commands.json .

./build/fluid --p-type=FLOAT --v-type=FAST_FIXED\(13,7\) --v-flow-type=FIXED\(32,5\)
