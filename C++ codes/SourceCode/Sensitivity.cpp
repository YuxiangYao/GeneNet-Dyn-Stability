#include <iostream>
#include <cstdio>
#include "Boolean_Function_Sensitivity.h"
using namespace std;
int main(int argc, char* argv[]){
    exe_B_Sensitivity_Analysis(argv);
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ Sensitivity.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o sensi.exe
    
SHELL: 
    #!/bin/bash
    ./sensi.exe R 1000 73001 & 		# Random BF
    ./sensi.exe E 1000 73002 & 		# Effective BF
    ./sensi.exe D 1000 73003 0 &	# Dominant BF
    ./sensi.exe D 1000 73004 9 &	# Dominant BF
    ./sensi.exe C 1000 73005 1 0 & 	# 0-Canalized BF
    ./sensi.exe C 1000 73006 1 9 & 	# Random Canalized BF
    ./sensi.exe C 1000 73007 2 0 & 	# 0-type 2-Nested Canalized BF
    ./sensi.exe C 1000 73008 2 9 & 	# Random 2-Nested Canalized BF
    ./sensi.exe P 1000 73009 0 &	# Post BF (only a^2)
    ./sensi.exe P 1000 73010 9 &	# Post BF (A^2+a^2)
    ./sensi.exe U 1000 73011 0 &	# monotone decreasing BF
    ./sensi.exe U 1000 73012 9 &	# Any monotone BF
    ./sensi.exe T 1000 73013 0 &	# 0-type Threshold BF
    ./sensi.exe T 1000 73014 9 &	# Random Threshold BF

FIGURE:
	Fig. 10a, 10b, S3.
*/ 