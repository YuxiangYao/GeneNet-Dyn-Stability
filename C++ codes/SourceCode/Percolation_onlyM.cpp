#include <iostream>
#include <cstdio> 
#include "Boolean_Model_MixedNK.h"
using namespace std;
int main(int argc, char* argv[]){
    exe_B_Percoaltion_Sample(argv);
    return 0;
}
/*
COMPILE: [...] is optional.
    g++ Percolation_onlyM.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o perco_m.exe
    
SHELL: 
    #!/bin/bash
    ./perco_m.exe U 4 2500 1000 25001 1 &	# monotone increasing BF & initial N0=N1
    ./perco_m.exe U 4 2500 1000 25002 0 &	# monotone decreasing BF & initial N0=N1

FIGURE: 
    Fig. 9.
*/
