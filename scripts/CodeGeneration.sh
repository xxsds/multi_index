#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: ./CodeGeneration n k"
    echo "Produces output files perm_n_k.hpp and perm_n_k.cpp"
      exit 1 
else
    python Permutations.py $1 $2 2> out && python CodeGeneration.py out perm_$1_$2
    mv perm_$1_$2.hpp ../include/multi_idx/perm_$1_$2.hpp
    mv perm_$1_$2.cpp ../lib/perm_$1_$2.cpp
    exit 0
fi
