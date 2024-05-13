#ifndef BOOLEAN_MODEL_MIXEDNK_H_INCLUDED
#define BOOLEAN_MODEL_MIXEDNK_H_INCLUDED
#include "Boolean_Function.h"
#include "Network_GraphTheory.h"
// class of mixed NK model analysis.
class B_Mixed_NK {
    public:
        char *ProjectName;
        int K;// K inputs
        int N;// System size
        int *TotalTruthTable;// Total truth tables (Length: 2^K*total)
        int *Parents;// Length: N*K
        int Loop;// Repeat times for analyzing
        int *ss;// System's s1,s2 old/new states
        int *Labels;// Label: all zeros to control
        double PQ[2];// Save bias PQ[0]:1, PQ[1]=1-PQ[1]
        int orders[15]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384};// [WARNING] K<15 !!!
        FILE *fpss;
        void Initialization(char **argv,char *FileName);
        void Delete();
        void SystemEvolution(int *sos,int Steps);// System evolution N-agent. [NOTE] ss is pionter_0 or _2N
        void B_Mixed_NK_MapSet(char **argv,double ratio);// Set system's N mapping tables. 
    protected:
        void B_Mixed_NK_R();
        void B_Mixed_NK_C(double ratios,char *CanaType,char *Random);
        void B_Mixed_NK_P(double ratios,char *PostType);
        void B_Mixed_NK_D(double ratios,char *DominType);
        void B_Mixed_NK_E(double ratios);
        void B_Mixed_NK_U(double ratios,char *UnateType);
        void B_Mixed_NK_T(double ratios,char *ThresType);
};
// Analysis of Derrida Curve.
class B_Mixed_NK_Derrida: public B_Mixed_NK {
    public:
        void InitialTwoVector(double distance,double bias);// Set two vectors with pointed distance.
        double TwoVectorDist();// Return two vectors normal distance.
        void SetParents_Derrida();// Set Derrida type parents.
};
// Analysis of stable state Percolation.
class B_Mixed_NK_Percolation: public B_Mixed_NK {
    public:
        int Length;
        int Height;
        NetworkSys tmp_net;
        void InitialOneVector();// Initial a ternary state.
        void SysEvoWindows(int Windows);// Set a window to observe "stable" nodes.
        void SetParents_Lattice();// Default lattice Para is K, length is sqrt(N);
        double LargestPercolationCluster();// Return largest cluster of percolation classes.
        double FractionStable();// Return fraction of stable nodes.
    protected:
        void TemporaryDirectedNet(struct Regulation *node);// Set a temporary directed network according to site's occupying.
        double AddressBottom2Top();// Address a path from bottom to top.
};

// Basic operations.
void B_Mixed_NK::Initialization(char **argv,char *FileName){// char *Type,char *Ks,char *SysSize,char *Times,char *RandSeed,
    K=atoi(argv[2]);// Ks
    N=atoi(argv[3]);// SysSize
    Loop=atoi(argv[4]);// Times
    mt.seed(atoi(argv[5]));// RandSeed
    TotalTruthTable=(int *)malloc(orders[K]*N*sizeof(int));// Total truth tables (Length: 2^K*total).
    ss=(int *)malloc(4*N*sizeof(int));// Quadruple, old.s1/s2 and new.s1/s2 
    Parents=(int *)malloc(K*N*sizeof(int));
    Labels=(int *)malloc(N*sizeof(int));//memset(Labels,0,N*sizeof(int));
    int ii;
    for(ii=0;ii<N;ii++)Labels[ii]=ii;// Set all node's code for mapping configuration.
    PQ[0]=PQ[1]=0.50;// {p,q} bias for random function.
    ProjectName=(char *)malloc(100*sizeof(char));
    char *cc[8],tmp[2]={"_"};strcpy(ProjectName,FileName);
    cc[0]=ProjectName;memcpy(cc+1,argv+1,5*sizeof(char *));
    StringPaste(cc,tmp,6);
    switch(argv[1][0]){
        case 'R':case 'E':ii=0;break;
        case 'D':case 'P':case 'U':case 'T':cc[1]=argv[6];ii=2;break;
        case 'C':cc[1]=argv[6];cc[2]=argv[7];ii=3;break;
        default:printf("Error in function type.\n");getchar();break;}
    StringPaste(cc,tmp,ii);
}
void B_Mixed_NK::Delete(){
    free(TotalTruthTable);free(ss);free(ProjectName);free(Labels);free(Parents);
    TotalTruthTable=ss=Labels=Parents=NULL;ProjectName=NULL;
}
void B_Mixed_NK::SystemEvolution(int *sos,int Steps){
    int ii,jj,kk,index,Length=orders[K],*old,*next,*changer,*maps,*pars;
    next=sos+N+N;old=sos;
    for(ii=0;ii<Steps;ii++){
        maps=TotalTruthTable;pars=Parents;
        for(jj=0;jj<N;jj++){
            index=0;
            for(kk=0;kk<K;kk++,pars++){// Default: all in-deg is same.
                index+=(*(old+*pars))*orders[kk];}
            next[jj]=maps[index];maps+=Length;}
        changer=old;old=next;next=changer;}
}
// Analyze various mixed NK-model.
void B_Mixed_NK::B_Mixed_NK_R(){
    int ii,*Maps=TotalTruthTable,Length;
    boolMap RBF;
    RBF.k=K;RBF.length=Length=orders[K];
    RBF.Initial();
    for(ii=0;ii<N;ii++,Maps+=Length){
        RBF.maptab=Maps;
        RBF.RBF_p(0.500);}
    RBF.Reset();
}
void B_Mixed_NK::B_Mixed_NK_E(double ratios){// ratios denotes specific TF class.
    int ii,part=(int)(ratios*N),Length;
    std::shuffle(Labels,Labels+N,mt);// Shuffle the node's code.
    boolMap_Effec EBF;
    EBF.k=K;EBF.length=Length=orders[K];
    EBF.Initial();
    for(ii=0;ii<part;ii++){// Specific TF part.
        EBF.maptab=TotalTruthTable+Length*Labels[ii];
        EBF.Effec_Gen();}
    for(ii=part;ii<N;ii++){
        EBF.maptab=TotalTruthTable+Length*Labels[ii];
        EBF.RBF_p(0.500);}
    EBF.Reset();
}
void B_Mixed_NK::B_Mixed_NK_D(double ratios,char *DominType){
    int ii,part=(int)(ratios*N),Length;
    std::shuffle(Labels,Labels+N,mt);// Shuffle the node's code.
    boolMap_Domin DBF;
    DBF.k=K;DBF.length=Length=orders[K];
    DBF.Initial();DBF.Domin_Type=atoi(DominType);
    DBF.In_Enumerate();DBF.Domin_Pre_Major();
    for(ii=0;ii<part;ii++){
        DBF.maptab=TotalTruthTable+Length*Labels[ii];
        DBF.Domin_Gen();}
    for(ii=part;ii<N;ii++){
        DBF.maptab=TotalTruthTable+Length*Labels[ii];
        DBF.RBF_p(0.500);}
    DBF.Reset();
}
void B_Mixed_NK::B_Mixed_NK_P(double ratios,char *PostType){
    int ii,part=(int)(ratios*N),Length;
    std::shuffle(Labels,Labels+N,mt);// Shuffle the node's code.
    boolMap_Post PBF;
    PBF.k=K;PBF.length=Length=orders[K];
    PBF.Initial();PBF.Post_Type=atoi(PostType);
    memset(PBF.labels,0,Length*sizeof(int));// Post_Ge
    for(ii=0;ii<part;ii++){
        PBF.maptab=TotalTruthTable+Length*Labels[ii];
        PBF.Post_Gen();}
    for(ii=part;ii<N;ii++){
        PBF.maptab=TotalTruthTable+Length*Labels[ii];
        PBF.RBF_p(0.500);}
    PBF.Reset();
}
void B_Mixed_NK::B_Mixed_NK_U(double ratios,char *UnateType){
    int ii,part=(int)(ratios*N),Length;
    std::shuffle(Labels,Labels+N,mt);// Shuffle the node's code.
    boolMap_Unate UBF;
    UBF.k=K;UBF.length=Length=orders[K];
    UBF.Initial();UBF.Unate_Type=atoi(UnateType);
    memset(UBF.labels,0,Length*sizeof(int));// Unate_Gen
    for(ii=0;ii<part;ii++){
        UBF.maptab=TotalTruthTable+Length*Labels[ii];
        UBF.Unate_Gen();}
    for(ii=part;ii<N;ii++){
        UBF.maptab=TotalTruthTable+Length*Labels[ii];
        UBF.RBF_p(0.500);}
    UBF.Reset();
}
void B_Mixed_NK::B_Mixed_NK_T(double ratios,char *ThresType){
    int ii,part=(int)(ratios*N),Length;
    std::shuffle(Labels,Labels+N,mt);// Shuffle the node's code.
    boolMap_Thres TBF;
    TBF.k=K;TBF.length=Length=orders[K];
    TBF.Initial();TBF.In_Enumerate();
    TBF.Thres_Type=atoi(ThresType);
    memset(TBF.labels,0,Length*sizeof(int));// Unate_Gen
    for(ii=0;ii<part;ii++){
        TBF.maptab=TotalTruthTable+Length*Labels[ii];
        TBF.Thres_Gen(0.500);}
    for(ii=part;ii<N;ii++){
        TBF.maptab=TotalTruthTable+Length*Labels[ii];
        TBF.RBF_p(0.500);}
    TBF.Reset();
}
void B_Mixed_NK::B_Mixed_NK_C(double ratios,char *CanaType,char *Random){
    int ii,part=(int)(ratios*N),Length;
    std::shuffle(Labels,Labels+N,mt);// Shuffle the node's code.
    boolMap_Cana CBF;
    CBF.k=K;CBF.length=Length=orders[K];
    CBF.Cana_Initial();CBF.In_Enumerate();
    CBF.Deep=atoi(CanaType);CBF.Fixed=atoi(Random);
    for(ii=0;ii<part;ii++){
        CBF.maptab=TotalTruthTable+Length*Labels[ii];
        CBF.Cana_Gen();}
    for(ii=part;ii<N;ii++){
        CBF.maptab=TotalTruthTable+Length*Labels[ii];
        CBF.RBF_p(0.500);}
    CBF.Cana_Reset();
}
void B_Mixed_NK::B_Mixed_NK_MapSet(char **argv,double ratios){
    switch(argv[1][0]){
        case 'R':B_Mixed_NK::B_Mixed_NK_R();break;
        case 'E':B_Mixed_NK::B_Mixed_NK_E(ratios);break;
        case 'D':B_Mixed_NK::B_Mixed_NK_D(ratios,argv[6]);break;
        case 'P':B_Mixed_NK::B_Mixed_NK_P(ratios,argv[6]);break;
        case 'U':B_Mixed_NK::B_Mixed_NK_U(ratios,argv[6]);break;
        case 'T':B_Mixed_NK::B_Mixed_NK_T(ratios,argv[6]);break;
        case 'C':B_Mixed_NK::B_Mixed_NK_C(ratios,argv[6],argv[7]);break;
        default:printf("Invalid Type.\n");getchar();break;}
}

// Derrida Curve for all specific BFs. 
void B_Mixed_NK_Derrida::InitialTwoVector(double distance,double bias){// [NOTE] Here, we only consider Initial Dis=0.10.
    int ii,tmp,Dist=(int)(N*distance);
    int *part1=ss,*part2=ss+N,*Selected=part2+N,*Flags=Selected+N;
    for(ii=0;ii<N;ii++,part1++){
        *part1=r_binary(mt);}
        // *part1=runif(mt)<bias;}
    part1=ss;
    memset(Flags,0,N*sizeof(int));
    if(distance<0.5){// Small part.
        memcpy(part2,ss,sizeof(int)*N);
        Choose_NK(N,Dist,Selected,Flags,-666);// Last Para no use here.
        for(ii=0;ii<Dist;ii++){
            tmp=Selected[ii];
            part2[tmp]=(0==part1[tmp]);}}// Inverse bits.
    else {// Large part.
        for(ii=0;ii<N;ii++){// Get inverse states.
            part2[ii]=(0==part1[ii]);}// Inverse bits.
        Choose_NK(N,N-Dist,Selected,Flags,-666);
        for(ii=0;ii<N-Dist;ii++){
            tmp=Selected[ii];
            part2[tmp]=part1[tmp];}}
}
double B_Mixed_NK_Derrida::TwoVectorDist(){
    int ii,*ipt1=ss,*ipt2=ss+N,sum=0;
    for(ii=0;ii<N;ii++,ipt1++,ipt2++){
        sum+=(*ipt1!=*ipt2);}
    return (double)(sum)/(double)(N);
}
void B_Mixed_NK_Derrida::SetParents_Derrida(){
    int ii,*pars=Parents;
    memset(ss,0,N*sizeof(int));// [Note] Borrow slot "ss" to set Selected, candiate.
    for(ii=0;ii<N;ii++,pars+=K){
        Choose_NK(N,K,pars,ss,ii);}
}
void exe_B_MixedNK_Derrida(char **argv){
    double frac=0,io[2],*distance;
    int ii=0,jj;char projects[20]={"b_MixedNK"};
    B_Mixed_NK_Derrida DerridaCurve;
    DerridaCurve.Initialization(argv,projects);
    distance=(double*)malloc(DerridaCurve.Loop*sizeof(double));
    DerridaCurve.fpss=fopen(strcat(DerridaCurve.ProjectName,".txt"),"w");
    for(ii=0;ii<101;ii++){
        for(jj=0;jj<DerridaCurve.Loop;jj++){
            DerridaCurve.SetParents_Derrida();
            DerridaCurve.B_Mixed_NK_MapSet(argv,frac);
            DerridaCurve.InitialTwoVector(0.10,0.80);// Distance keeps 0.1000, 0.80 only takes effect in removing comment Line 214.
            DerridaCurve.SystemEvolution(DerridaCurve.ss,1000);// S1 vector
            DerridaCurve.SystemEvolution(DerridaCurve.ss+DerridaCurve.N,1000);// S2 vector
            distance[jj]=DerridaCurve.TwoVectorDist();}
        Mean_SD(distance,DerridaCurve.Loop,io);
        fprintf(DerridaCurve.fpss,"%.2f,%.4f,%.4f\n",frac,io[0],io[1]);
        frac+=0.01;}
    fclose(DerridaCurve.fpss);free(distance);
    DerridaCurve.Delete();
}
// Derrida Curve for only dominating class.
void exe_B_MixedNK_D_verify(char **argv){
    double frac=0.80,*distance;// Only set 0.8 fracton, d_0=0.1, N=1000.
    int jj;char projects[20]={"b_MixedNK_D_distr"};
    B_Mixed_NK_Derrida DerridaCurve;
    DerridaCurve.Initialization(argv,projects);
    distance=(double*)malloc(DerridaCurve.Loop*sizeof(double));
    DerridaCurve.fpss=fopen(strcat(DerridaCurve.ProjectName,".txt"),"w");
    for(jj=0;jj<DerridaCurve.Loop;jj++){
        DerridaCurve.SetParents_Derrida();
        DerridaCurve.B_Mixed_NK_MapSet(argv,frac);
        DerridaCurve.InitialTwoVector(0.10,0.80);// Distance keeps 0.1000
        DerridaCurve.SystemEvolution(DerridaCurve.ss,1000);// S1 vector
        DerridaCurve.SystemEvolution(DerridaCurve.ss+DerridaCurve.N,1000);// S2 vector
        distance[jj]=DerridaCurve.TwoVectorDist();
        fprintf(DerridaCurve.fpss,"%.4f\n",distance[jj]);}
    fclose(DerridaCurve.fpss);free(distance);
    DerridaCurve.Delete();
}
// Percolation lattice.
void B_Mixed_NK_Percolation::InitialOneVector(){
    for(int ii=0;ii<N;ii++){
        ss[ii]=r_binary(mt);}
        //ss[ii]=runif(mt)<0.8;}
}
void B_Mixed_NK_Percolation::SysEvoWindows(int Windows){
    int ii,jj,kk,index,*old,*next,*changer,*maps,*pars,*jumper,*single_s;
    int Length=orders[K];
    old=ss;next=ss+N;jumper=next+N;single_s=jumper+N;
    memset(jumper,0,N*sizeof(int));
    for(ii=0;ii<Windows;ii++){
        maps=TotalTruthTable;pars=Parents;
        for(jj=0;jj<N;jj++){
            index=0;
            for(kk=0;kk<K;kk++,pars++){// Default: all in-deg is same.
                index+=(*(old+*pars))*orders[kk];}
            next[jj]=maps[index];maps+=Length;
            jumper[jj]+=(old[jj]!=next[jj]);}// Count each node's changing in observing window.
        changer=old;old=next;next=changer;}
    for(ii=0;ii<N;ii++){// If the node keep stable in the window?
        single_s[ii]=(jumper[ii]==0);}
}
void B_Mixed_NK_Percolation::SetParents_Lattice(){
    int ii,*pars=Parents;
    for(ii=0;ii<N;ii++,pars+=K){
        Auxiliary_LT_U(Length,Height,K,ii,pars);}
}
void B_Mixed_NK_Percolation::TemporaryDirectedNet(struct Regulation *node){
    int *if_stable=ss+N+N+N;
    int id=node->code;// This network is unidirectional.
    int ups,downs,lefts,rights;
    ups=id+Length;downs=id-Length;
    if(id%Length==0){lefts=id+Length-1;}else {lefts=id-1;}// Left periodic boundary condition
    if(id%Length==Length-1){rights=id-Length+1;}else {rights=id+1;}// Right periodic boundary condition
    if(downs>=0&&if_stable[downs]>0)InsertChildsNode(node,downs);// Down has boundary
    if(ups<N&&if_stable[ups]>0)InsertChildsNode(node,ups);// Up has boundary 
    if(if_stable[lefts]>0){InsertChildsNode(node,lefts);}
    if(if_stable[rights]>0){InsertChildsNode(node,rights);}
}
double B_Mixed_NK_Percolation::AddressBottom2Top(){
    int ii,logi=0;
    double results;
    int *Labels=ss+N;
    int *routers=tmp_net.InDeg+3*N;// Here is routetable.
    int *ipt=routers+N-Length;
    for(ii=0;ii<Length;ii++){
        if(*(ipt+ii)>N);// Has a path? (Only check top line)
        else {logi=1;break;}}
    if(logi){// Find all nodes in this path.
        ipt=routers;results=0;
        for(ii=0;ii<N;ii++){
            results+=*(ipt++)<N;}}// Default: 2*N
    else results=0;
    if(results>0){// Remove some common largest candidate labels.
        for(ii=0;ii<Length;ii++){
            if(routers[ii]<N){
                Labels[ii]=0;}}}
    return results;
}
double B_Mixed_NK_Percolation::LargestPercolationCluster(){
    int ii; double sum,tmp;
    int *if_stable=ss+N+N+N;// [Note] After window step, ss+2N is the slot.
    //int *routetables=ss;// [Note] Borrow ss slot, ss+0.
    int *ReduceLabel=ss+N;// [Note] Borrow ss slot, ss+N.
    for(ii=0;ii<N;ii++){// Set a temporary directed network (Only left-right periodic).
        if(if_stable[ii]>0)B_Mixed_NK_Percolation::TemporaryDirectedNet(tmp_net.Network+ii);}
    sum=0;
    for(ii=0;ii<Length;ii++)ReduceLabel[ii]=if_stable[ii];
    for(ii=0;ii<Length;ii++){
        if(ReduceLabel[ii]>0){
            tmp_net.ShortestPath(ii);// The routetables saved in [tmp_net.InDeg+N*3] location.
            tmp=B_Mixed_NK_Percolation::AddressBottom2Top();// Find the pathways.
            if(tmp>sum)sum=tmp;}}
    tmp_net.Reset_Network();// Reset the network framework.
    return sum;
}
double B_Mixed_NK_Percolation::FractionStable(){
    int ii,sum=0,*nodes=ss+N+N+N;
    for(ii=0;ii<N;ii++){
        sum+=nodes[ii]>0;}
    return ((double)sum)/((double)N);
}
void exe_B_MixedNK_Percolation(char **argv){
    double frac=0,io[2],*perco,*stabler,Perco_Ratio; 
    int ii=0,jj,n_loop,n_size,*sss;char projects[20]={"b_Perco"};
    B_Mixed_NK_Percolation StablePerco;
    StablePerco.Initialization(argv,projects);
    StablePerco.Height=StablePerco.Length=(int)(sqrt(StablePerco.N));
    perco=(double*)malloc(StablePerco.Loop*sizeof(double));
    stabler=(double*)malloc(StablePerco.Loop*sizeof(double));
    StablePerco.tmp_net.total=StablePerco.N;
    StablePerco.tmp_net.Initial_Network();// Initial network.
    sss=StablePerco.ss;n_loop=StablePerco.Loop;n_size=StablePerco.N;
    StablePerco.fpss=fopen(strcat(StablePerco.ProjectName,".txt"),"w");
    for(ii=0;ii<101;ii++){
        Perco_Ratio=0;
        for(jj=0;jj<StablePerco.Loop;jj++){
            StablePerco.SetParents_Lattice();
            StablePerco.B_Mixed_NK_MapSet(argv,frac);
            StablePerco.InitialOneVector();// A random ternary vector.
            StablePerco.SystemEvolution(sss,1000);// Only S1 vector simulation
            StablePerco.SysEvoWindows(500);// Observing window [Ref. J. theor. Biol. (1988) 135, 255-261]
            perco[jj]=StablePerco.LargestPercolationCluster()/n_size;
            stabler[jj]=StablePerco.FractionStable();
            Perco_Ratio+=perco[jj]>0;}
        Mean_SD(perco,StablePerco.Loop,io);
        fprintf(StablePerco.fpss,"%.2f,%.4f,%.4f,%.4f,",frac,Perco_Ratio/n_loop,io[0],io[1]);
        Mean_SD(stabler,StablePerco.Loop,io);
        fprintf(StablePerco.fpss,"%.4f,%.4f\n",io[0],io[1]);
        frac+=0.01;}
    StablePerco.tmp_net.Delete_Network();
    fclose(StablePerco.fpss);
    StablePerco.Delete();
}
// Show one sample of percolation of U-type hybrid models.
void exe_B_Percoaltion_Sample(char **argv){
    double frac=1.000,perco,stabler,Perco_Ratio; 
    int ii=0,jj,n_loop,n_size,*sss,*ttt;
    char projects[20]={"b_Perco_Sample"};
    B_Mixed_NK_Percolation StablePerco;
    StablePerco.Initialization(argv,projects);
    StablePerco.Height=StablePerco.Length=(int)(sqrt(StablePerco.N));
    StablePerco.tmp_net.total=StablePerco.N;
    StablePerco.tmp_net.Initial_Network();// Initial network.
    sss=StablePerco.ss;ttt=StablePerco.ss+StablePerco.N*3;
    n_loop=StablePerco.Loop;n_size=StablePerco.N;
    StablePerco.fpss=fopen(strcat(StablePerco.ProjectName,".txt"),"w");
    for(jj=0;jj<StablePerco.Loop;jj++){
        StablePerco.SetParents_Lattice();
        StablePerco.B_Mixed_NK_MapSet(argv,frac);
        StablePerco.InitialOneVector();// A random ternary vector.
        StablePerco.SystemEvolution(sss,1000);// Only S1 vector simulation
        StablePerco.SysEvoWindows(500);// Observing window [Ref. J. theor. Biol. (1988) 135, 255-261]
        perco=StablePerco.LargestPercolationCluster()/n_size;
        if(perco>0){// Occur percolation
            for(ii=0;ii<StablePerco.N;ii++){
                if(ttt[ii]){// Stable
                    fprintf(StablePerco.fpss,"%d\n",sss[ii]);}
                else {// -1 means not stable, to distinguish from stable-0, stable-1
                    fprintf(StablePerco.fpss,"-1\n");}}
            jj=100000;}// If finding a percolation, stop simulations.
        else ;}
    StablePerco.tmp_net.Delete_Network();
    fclose(StablePerco.fpss);
    StablePerco.Delete();
}
#endif // BOOLEAN_MODEL_MIXEDNK_H_INCLUDED
