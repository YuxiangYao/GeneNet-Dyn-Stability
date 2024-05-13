#include <iostream>
#include <cstdio>
#include <malloc.h> 
#include "Random_Setting.h"
#include "Analysis_BioNets_Linux.h"
//#include "Analysis_BioNets_Windows.h"
using namespace std;
int main(int argc, char* argv[]){
    // Analyze sensitivity of threshold based models.
    char dir[]={"//change//as//your//dir//"},lists[]={"list_L.txt"};
    char res1[]={"LBN_AvSen.txt"}, res2[]={"LBN_FunFrac.txt"};
    exe_AvSensitivity_Logic(dir,lists,res1);
    // Identify the logical-based Nets can be also which types of BF.
    exe_FunFracComp_Logic(dir,lists,res2);
    return 0;
}
/*  
COMPILE: [...] is optional.
    g++ Analysis_L_BioNet.cpp [-L/your/z3lib/dir/bin] [-I /your/z3lib/dir/include] -lz3 -o bionet_BF.exe
    ./bionet_BF.exe

FIGURE:
    Fig. 2a, 2c, 2d, 3a.

MODIFY:
    Change the header file according to the system you employ. [Line5 <=> Line6]

NOTE: 
    Network list file (list_L.txt), and series Net-folders should be in the same folder/path.
*/
