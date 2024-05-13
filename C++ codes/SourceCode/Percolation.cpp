#include <iostream>
#include <cstdio> 
#include "Boolean_Model_MixedNK.h"
using namespace std;
int main(int argc, char* argv[]){
	exe_B_MixedNK_Percolation(argv);
    return 0;
}
/*
COMPILE: [...] is optional.
    g++ Percolation.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o perco.exe
    
SHELL: 
    #!/bin/bash
    ./perco.exe R 1000 77001 & 		# Random BF
    ./perco.exe E 1000 77002 & 		# Effective BF
    ./perco.exe D 1000 77003 0 &	# Dominant BF
    ./perco.exe D 1000 77004 9 &	# Dominant BF
    ./perco.exe C 1000 77005 1 0 & 	# 0-Canalized BF
    ./perco.exe C 1000 77006 1 9 & 	# Random Canalized BF
    ./perco.exe C 1000 77007 2 0 & 	# 0-type 2-Nested Canalized BF
    ./perco.exe C 1000 77008 2 9 & 	# Random 2-Nested Canalized BF
    ./perco.exe P 1000 77009 0 &	# Post BF (only a^2)
    ./perco.exe P 1000 77010 9 &	# Post BF (A^2+a^2)
    ./perco.exe T 1000 77011 0 &	# 0-type Threshold BF
    ./perco.exe T 1000 77012 9 &	# Random Threshold BF
    ./perco.exe U 1000 77013 0 &	# monotone decreasing BF
    ./perco.exe U 1000 77014 9 &	# Any monotone BF
    ./perco.exe U 1000 77015 1 &	# monotone increasing BF

FIGURE: 
    Fig. 8.

MODIFY:
    For the last program, before compile, please comment out Line 286 and recover Line 287.
    ./perco.exe D 1000 77016 9 &	# Dominant BF, Initial bias of "state" is 0.80.
*/
