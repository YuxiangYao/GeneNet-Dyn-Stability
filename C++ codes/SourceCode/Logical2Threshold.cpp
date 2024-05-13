#include <iostream>
#include <cstdio>
#include <malloc.h> 
#include "Random_Setting.h"
#include "Analysis_BioNets_Linux.h"
//#include "Analysis_BioNets_Windows.h"
using namespace std;
int main(int argc, char* argv[]){
    // Test: convert logical to threshold.
    char dir[]={"//change//as//your//dir//"},lists[]={"list_L.txt"};
    char res0[]={"Logic2Thres.txt"};
    exe_Logic_to_Threshold(dir,lists,res0);
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ Analysis_L_BioNet.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o Logic2Thres.exe
    ./Logic2Thres.exe

FIGURE:
    None.

MODIFY:
    Change the header file according to the system you employ. [Line5 <=> Line6]

NOTE: 
    Network list file (list_L.txt), and series Net-folders should be in the same folder/path.
*/
