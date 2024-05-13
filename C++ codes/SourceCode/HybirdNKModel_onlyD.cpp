#include <iostream>
#include <cstdio> 
#include "Boolean_Model_MixedNK.h"
using namespace std;
int main(int argc, char* argv[]){
    exe_B_MixedNK_D_verify(argv);// for Fig.6b
    //exe_B_MixedNK_Derrida(argv);// for Fig.6c
    return 0;
}
/*
MODIFY:
    Before complie, comment out Line 213 and remove the comment of Line 214 for Fig.6b
    Comment out Line 6 in this .cpp and recover Line 7 for Fig.6c. (Also need Line 213 -> Line 214)

COMPILE: [...] is optional.
    g++ HybirdNKModel_D.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o derr_d.exe
    
SHELL-1: Fig.6b
    #!/bin/bash
    ./derr_d.exe D 1000 67321 9 &		# Dominant BF, Initial Dist=0.5 (Change the value in the header file Line 275)
    ./derr_d.exe D 1000 89032 9 &		# Dominant BF, Initial Dist=0.9 (Change the value in the header file Line 275)

SHELL-2: Fig.6c
    #!/bin/bash
    ./derr_d.exe D 1000 89331 1 &		# Dominant BF
    ./derr_d.exe D 1000 89332 9 &		# Dominant BF

FIGURE: 
    Fig. 6b,6c
*/