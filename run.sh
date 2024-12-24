#!/bin/bash

mkdir -p build 

cmake -S . \
      -B build \
      -G Ninja \
      -DCMAKE_CXX_COMPILER=clang++-18 \
      -DCMAKE_C_COMPILER=clang-18 \
      -DTYPES="FAST_FIXED(64,8)" \
      -DSIZES="S(36,84)" \
      -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping" \
      -DCMAKE_BUILD_TYPE=release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=1 
#-DSIZES="S(36,84),S(10,10),S(42,1337)"
#-DTYPES="FLOAT,FAST_FIXED(13,7),FAST_FIXED(32,5),DOUBLE"
#-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -Wall -Wextra -pedantic -fprofile-instr-generate -fcoverage-mapping -O2" \
cmake --build build --parallel $(nproc) || exit 1
cp -f build/compile_commands.json .

#./build/fluid --p-type=FLOAT --v-type=FAST_FIXED\(32,5\) --v-flow-type=FAST_FIXED\(32,5\) --field="fields/example1"
./build/fluid --p-type="FAST_FIXED(64,8)" --v-type="FAST_FIXED(64,8)" --v-flow-type=FAST_FIXED"(64,8)" --field="fields/example1" --j=2
#./build/fluid-stupid