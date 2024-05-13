# GeneNet-Dyn-Stability
## **Introduction**
C++ code and data files of "Canalization and competition: the cornerstone of genetic network’s dynamic stability and evolution". The detail descriptions of algorithm, setting, and analysis can be found in the paper. 

All simulations are implemented by these C++ codes ( $\geqslant$ `CPP98`) and depicted by `R` (v4.2.3) with `ggplot2` (v3.4.4). 

## **Code & File**
### **Simulation and analysis codes**
There are two folders, HeaderFile and SourceCode. `HeaderFile` contains all analysis functions' declarations and definitions (All are `.h` files; No individual `.cpp` files for sake of modularity). `SourceCode` contains all source codes to compile (Only `.cpp` files), which can output the executable files to conduct analysis and simulation. The detail parameter configurations, annotations, and pseudo-shell codes are also recorded in corresponding `.cpp` files. Compling one `.cpp` file needs to ensure corresponding header files and the source code in same paths. All file's dependencies are appropriate. Brief introduction of source codes as follow,

-------
- **Sensitivity.cpp** :  Analyze the sensitivity of various kinds of speical Boolean function (BF). 
- **BF_energy.cpp** :  Return the total/relative energy of BF.
-------
- **HybirdNKModel.cpp** : Simulate the hybird NK model to obtain calculate the Derrida Curves.
- **HybirdNKModel_onlyD.cpp** : Simulate the hybird NK model for "case Dominant BF" with special settings.
- **Percolation.cpp** : Observe percolation of various hybird NK model.
- **Percolation_onlyM.cpp** : Observe percolation of hybird NK model for "case Monotone BF" with special settings.
-------
- **Analysis_L_BioNet.cpp** : Analyze bionet's properties and SBF's categories and of logical-based networks.
- **Analysis_T_BioNet.cpp** : Analyze bionet's properties and SBF's categories and of threshold-based networks.
- **Logical2Threshold.cpp** ：Check the logical-based BFs converted to threshold-based ones.
- **BioNetDyn.cpp** : Coarse-grained analyze bionet's attractor by dynamic simulation.
-------
- **Evolution.cpp** : Simulate the evolutionary process and reveal the C/T BF's capacity.
- **Evolution_GenPoint.cpp** : Only analyze the attractor's cases at the 50/100/150 generation.


### **Genetic network files**
- **LogNet.zip**: Includes 71 logical-based genetic networks from Cell Collective Website. Each folder represents one specific network. `xxx.csv` file is the mapping table of gene `xxx`. Mapping table formats instead of logic expression ones is employed for the sake of facilitating downstream analysis.
- **ThrNet.zip**: Includes 8 threshold-based genetic networks from published literatures. Unlike logical-based nets, each threshold-based net is only consist of two paired files, `xxx.ids` and `xxx.topo`. `xxx.ids` records the genenames and identifier. `xxx.topo` save the topological and regulatory relations (Third column: 1 for activate, 2 for inhibit). 

> [!NOTE]
> Some codes utilize a third-party library (related to framework switching parts), Z3solver (https://github.com/Z3Prover/z3/releases). Please ensure correct installation and add the path in the `LD_LIBRARY_PATH`. Command `g++ your.cpp [...] -lz3 -o yours.exe` to compile program.

> [!TIP]
> Some simulations are time-consuming (Such as attractor finding or percolation). You can use `OpenMP` or `pthread` multi-processing methods or split programs by shells to speed up the numerical analysis. 

## **Questions** 
For any question about the program, please contact Dr. Yuxiang Yao, Email: y.x.yao@foxmail.com
