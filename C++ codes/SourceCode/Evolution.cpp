#include <iostream>
#include <cstdio>
#include <malloc.h> 
#include <string.h>
#include "Random_Setting.h"
#include "Evolution.h"
using namespace std;
int main(int argc, char* argv[]){
    // Argu: 1:FileName, 2:RepeatTimes, 3:RandSeed
    char dirs[]={"//change//as//your//dir//"};
    // Here freely evolve the network.
    // Each probability of Swap-Edge, Map Reset, Add-Node, Del-Node, Add-Edge, Del-Edge
    double probs[]={1.0/6.0,1.0/2.0,0.5,2.0/3.0,5.0/6.0};// cumulative probability
    mt.seed(atoi(argv[3]));
    Evoluting_OnlyEdge_Map(dirs,argv[1],atoi(argv[2]),probs,1);
    
    // Here only expand network by adding node/edge 
    //mt.seed(atoi(argv[3]));
    //double probs[]={0.25,0.50,0.75,0.75,1.0};// cumulative probability
    //Evoluting_OnlyEdge_Map(dirs,argv[1],atoi(argv[2]),probs,1);
    
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ Evolution.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o Evolution.exe

FIGURE:
    Fig. 3b, 3c.

SHELL: 
    #!/bin/bash
    ./Evolution.exe evo001 1000 12345 & # name:evo001, 1000 repeats, random seed: 12345
    .....

MODIFY:
    * Change the header file according to the system you employ. [BioNets_Dynamic.h: Line3 <=> Line4]
    * If observing expanding network, please change the comments [this file: Line13~15 <=> Line18~20]

NOTE: 
    * Switch the mode: only C/T, C/T+random, only random, please change the Line695 in Evolution.h
*/