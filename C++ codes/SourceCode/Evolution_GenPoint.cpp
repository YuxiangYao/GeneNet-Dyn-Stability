#include <iostream>
#include <cstdio>
#include <malloc.h> 
#include <string.h>
#include "Random_Setting.h"
#include "Evolution.h"
using namespace std;
int main(int argc, char* argv[]){
    // Argu: 1:FileName, 2:Generation, 3:RandSeed
    char dirs[]={"//change//as//your//dir//"};
    
    // Here freely evolve the network.
    // Each probability of Swap-Edge, Map Reset, Add-Node, Del-Node, Add-Edge, Del-Edge
    double probs[]={1.0/6.0,1.0/2.0,0.5,2.0/3.0,5.0/6.0};// cumulative probability
    mt.seed(atoi(argv[3]));
    Evoluting_EvoAttrLoop(dirs,argv[1],probs,atoi(argv[2]));
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ Evolution.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o Evolution.exe

FIGURE:
    Fig. 5a~5c.

SHELL: 
    #!/bin/bash
    ./Evolution.exe evo050_000 50 11001 &
    ./Evolution.exe evo050_002 50 11002 &
                .....
    ./Evolution.exe evo050_999 50 11999 &
    ./Evolution.exe evo050_999 50 11999 &
    ./Evolution.exe evo100_999 100 12000 &
                .....

MODIFY:
    * Change the header file according to the system you employ. [BioNets_Dynamic.h: Line3 <=> Line4]

NOTE: 
    * Switch the mode: only C/T, C/T+random, only random, please change the Line 695 in Evolution.h
    * SHELL employs parallel methods to speed up.
*/