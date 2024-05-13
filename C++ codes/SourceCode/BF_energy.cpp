#include <iostream>
#include <cstdio>
#include "Boolean_Function_Energy.h"
using namespace std;
int main(int argc, char* argv[]){
    exe_B_Energy_Analysis(argv);
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ BF_energy.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o energy.exe
    
SHELL: 
    #!/bin/bash
    ./energy.exe D 1000 75003 0 &	    # Dominant BF
    ./energy.exe D 1000 75004 9 &	    # Dominant BF
    ./energy.exe C 1000 75005 1 0 & 	# 0-Canalized BF
    ./energy.exe C 1000 75006 1 9 & 	# Random Canalized BF
    ./energy.exe C 1000 75007 2 0 & 	# 0-type 2-Nested Canalized BF
    ./energy.exe C 1000 75008 2 9 & 	# Random 2-Nested Canalized BF
    ./energy.exe P 1000 75009 0 &	    # Post BF (only a^2)
    ./energy.exe P 1000 75010 9 &	    # Post BF (A^2+a^2)
    ./energy.exe U 1000 75011 0 &	    # monotone decreasing BF
    ./energy.exe U 1000 75012 9 &	    # Any monotone BF
    ./energy.exe T 1000 75013 0 &	    # 0-type Threshold BF
    ./energy.exe T 1000 75014 9 &	    # Random Threshold BF

FIGURE:
    Fig. 10c,S4
*/ 
