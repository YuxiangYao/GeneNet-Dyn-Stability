#include <iostream>
#include <cstdio> 
#include "Boolean_Model_MixedNK.h"
using namespace std;
int main(int argc, char* argv[]){
    exe_B_MixedNK_Derrida(argv);
    return 0;
}
/*
COMPILE: [...] is optional.
    g++ HybirdNKModel.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o derr.exe
    
SHELL: (~75min)
    #!/bin/bash
    ./derr.exe R 1000 73001 & 		# Random BF
    ./derr.exe E 1000 73002 & 		# Effective BF
    ./derr.exe D 1000 73003 0 &		# Dominant BF
    ./derr.exe D 1000 73004 9 &		# Dominant BF
    ./derr.exe C 1000 73005 1 0 & 	# 0-Canalized BF
    ./derr.exe C 1000 73006 1 9 & 	# Random Canalized BF
    ./derr.exe C 1000 73007 2 0 & 	# 0-type 2-Nested Canalized BF
    ./derr.exe C 1000 73008 2 9 & 	# Random 2-Nested Canalized BF
    ./derr.exe P 1000 73009 0 &		# Post BF (only a^2)
    ./derr.exe P 1000 73010 9 &		# Post BF (A^2+a^2)
    ./derr.exe U 1000 73011 0 &		# monotone decreasing BF
    ./derr.exe U 1000 73012 9 &		# Any monotone BF
    ./derr.exe T 1000 73013 0 &		# 0-type Threshold BF
    ./derr.exe T 1000 73014 9 &		# Random Threshold BF

FIGURE: 
    Fig. 6a.
*/