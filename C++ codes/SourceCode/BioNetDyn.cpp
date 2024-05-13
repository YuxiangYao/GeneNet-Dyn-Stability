#include <iostream>
#include <cstdio>
#include <malloc.h> 
#include "BioNets_Dynamic.h"
using namespace std;
int main(int argc, char* argv[]){
    char dir[]={"//change//as//your//bionet//dir//"};
    // 1:Net ID, 2:OutDir, 3:RandSeed, 4-FileName, 5-RepaetTimes
    mt.seed(atoi(argv[3]));
    exe_SearchAttractor(dir,argv[1],argv[2],argv[4],atoi(argv[5]));
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ BioNetDyn.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o BioNetDyn.exe
    ./BioNetDyn.exe

FIGURE:
    Fig. 4b, 4c, 4d.

SHELL: 
    #!/bin/bash
    ./BioNetDyn.exe [NetID] [outdir] [Seed] [file index name] [1000] &
    .....
    (Employed seeds: 230624 -> 230694)

MODIFY:
    Change the header file according to the system you employ. [BioNets_Dynamic.h: Line3 <=> Line4]

NOTE: 
    * Before executing, make sure the folder [argv2:OutDir] avaiable.
    * Some sizes of BioNet are large, need more repeats and it is time-consuming. You can split one
      network's simulation by parallel methods.
*/