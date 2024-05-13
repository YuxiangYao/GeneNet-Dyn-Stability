/*
  Copyright (c) 2023-2024, YAO Yuxiang
    date: 2024-01-16
    version: 1.0
    
  [Use procedure] 
  [R] Random BF: Initial -> (In_Enumerate) -> { RBF_p } -> Reset
  [B] Block BF: Initial+Argu -> In_Enumerate -> Block_Pre_Major -> { Block_Gen } -> Reset
  [E] Effec BF: Initial -> (In_Enumerate) -> { Effec_Gen } -> Reset
  [C] Cana  BF: Cana_Initial+Argu -> In_Enumerate -> {Cana_Gen} -> Cana_Reset
  [P] Post  BF: Initial+Argu -> (In_Enumerate) -> {Post_Gen} -> Reset
  [T] Thres BF: Initial -> In_Enumerate -> Thres_Gen -> Reset
  [U] Unate BF: Initial+Argu -> (In_Enumerate) -> { Unate_Gen } -> Reset
    -> (*) means this step is not essential
    -> Must follow above orders to generate BFs and reset. If only repeating simulations rather than self-defined using, one can ignore above settings. 
    -> In order to clarify the function of each BF, we did not use overloading methods, thus each type BF inherit main Boolean function and define subclass.

  [References]
  Special classes of Bf can be found in:
  1) Canalized: The origins of order: self-organization
  2) Post: 10.1073/pnas.1534782100
  3) Block: in this paper
  4) Effective: 10.1006/jtbi.2002.3081
  5) Unate: 1398.10.1007/s11538-008-9304-7
  6) Threshold: 10.1103/PhysRevE.86.026114
*/
#ifndef BOOLEAN_FUNCTION_H_INCLUDED
#define BOOLEAN_FUNCTION_H_INCLUDED
#include "Common_Used.h"
#include <vector>
#include <string>
#include <z3++.h>

class boolMap {
    public:
        int k;// Number of input variables
        int length;// Length of mapping table
        int *maptab;// Mapping table (a vector, ordered by their input vecttors, like 000, 001, 010,...)
        int *enum_dict;// Mapping table dictionary
        int *labels;// For subclass generators usage
        // Functions & Operations
        void Initial();// Initialization.
        void Reset();// Reset.
        void RBF_p(double bias);// Return a Bernoulli distribution results with bias (p).
        void In_Enumerate();// Return a k-input enumeration dictionary (length: k*2^k).
        void Show_Map();// Show the truth table.
        double Num_1();// Count number of 1 in truth table.
        double Sensitivity();// Calculate sensitivity of BFs.
        double Energy();// Calculate function energy.
        void Exchange_Map();// Random exchange mapping table's values.
};
// Derived class: canalized BF
class boolMap_Cana: public boolMap {
    public:
        int Deep;// Canalized deep/level.
        int Fixed;// Should be fixed or random?
        int *Selected;// Which variables should be set?
        int *Canalzing;// In values.
        int *Canalized;// Out values.
        int is_Cana();// Judge mapping is Canalized class or not.
        void Cana_Initial();// Initialization
        void Cana_Reset();// Delete 
        void Cana_Gen();// General generator.
        void Cana_Pure_Gen();// Canalized class generator.
        void Cana_Nest_Gen();// Nest-canalized class generator (layer=2).
};
// Derived class: Post2 BF
class boolMap_Post: public boolMap {
    public: // [NOTE] Post class only consider {A2} or {a2} subclass.
        int Post_Type;
        int is_Post();// Judge mapping is Post class or not.
        void Post_Gen();// Post class generator.
};
// Derived class: dominant BF
class boolMap_Domin: public boolMap {
    public:
        int Domin_Type;
        int is_Domin();// Judge mapping is Block class or not.
        void Domin_Pre_Major();// Obtain enumerate dictionary major state.
        void Domin_Gen();// Block class generator.
};
// Derived class: effective BF
class boolMap_Effec: public boolMap {
    public:
        int is_Effec();// Judge mapping is effective class or not.
        void Effec_Gen();// Effective class generator.
};
// Derived class: monotone BF
class boolMap_Unate: public boolMap {
    public:
        int Unate_Type;// Set type of monotone BF, 1/0 is increase/descreas.
        int is_Unate();// Judge mapping is monotone class or not.
        int is_Mixed_Unate();// Judge mapping is monotone class or not.
        void Unate_Gen();// Monotone class generator.
};
// Derived class: threshold BF
class boolMap_Thres: public boolMap {
    public:
        int Thres_Type;
        int is_Thres();// Judge mapping is threshold class or not.
        void Thres_Gen(double p_bias);// Threshold class generator.
};

// Reset the instance.
void boolMap::Initial(){
    enum_dict=(int *)malloc(k*length*sizeof(int));
    labels=(int *)malloc(length*sizeof(int));
}
void boolMap::Reset(){
    k=length=0;free(labels);free(enum_dict);
    maptab=labels=enum_dict=NULL;
}
// Return a Bernoulli distribution results with bias (p).
void boolMap::RBF_p(double bias){
    for(int ii=0;ii<length;ii++)
        maptab[ii]=(runif(mt)<bias);
}
// Return a k-input enumeration dictionary (length: k*2^k).
void boolMap::In_Enumerate(){
    int ii,jj,index,bits[16]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384};
    for(ii=0;ii<length;ii++){
        index=ii;
        for(jj=0;jj<k;jj++,index+=length){
            enum_dict[index]=(ii&bits[jj])>0;}}
}
// Show the truth table.
void boolMap::Show_Map(){
    int ii,jj,index;
    for(ii=0;ii<length;ii++){
        std::cout << '[';
        index=ii+length*(k-1);
        for(jj=0;jj<k;jj++,index-=length){
            std::cout << enum_dict[index] << ' ';}
        std::cout << "] [" << maptab[ii] << ']' << std::endl;}
}
// Count number of 1 in truth table.
double boolMap::Num_1(){
    int ii;double num_1=0;
    for(ii=0;ii<length;ii++){
        num_1+=(maptab[ii]>0);}
    return num_1;
}
// Calculate sensitivity of BFs.
double boolMap::Sensitivity(){
    int ii,jj,kk,*ipt1,*ipt2,boundary_l=1,boundary_u=length>>1,sumer=0;
    for(ii=0;ii<k;ii++){
        ipt1=maptab;
        ipt2=maptab+boundary_l;
        for(jj=0;jj<boundary_u;jj++){
            for(kk=0;kk<boundary_l;kk++){
                sumer+=(ipt1[kk])!=(ipt2[kk]);}
            ipt1+=(boundary_l<<1);
            ipt2+=(boundary_l<<1);}
        boundary_l=boundary_l<<1;
        boundary_u=boundary_u>>1;}
    return (sumer<<1)/(double)(length);// [NOTE] Ergodic checking only half part, should be doubled.
}
// Calculate function energy.
double boolMap::Energy(){
    int ii,jj,code,energy=0,spin[length],bits[k];
    bits[0]=1;
    for(ii=1;ii<k;ii++)bits[ii]=bits[ii-1]<<1;
    for(ii=0;ii<length;ii++)spin[ii]=(maptab[ii]<<1)-1;
    for(ii=0;ii<length;ii++){
        for(jj=0;jj<k;jj++){
            code=ii^bits[jj];
            energy+=spin[ii]*spin[code];}}
    return (double) energy/2;
}
// Random exchange mapping table's values.
void boolMap::Exchange_Map(){
    int ii,id1,id2,exchanger,Less=length-1;
    for(ii=0;ii<(length<<1);ii++){// Should exchange "0.5*length" times.
        id1=(int)(runif(mt)*length);
        id2=(int)(runif(mt)*Less)+1;
        id2=(id1+id2)%length;// Get the loop.
        exchanger=maptab[id1];
        maptab[id1]=maptab[id2];
        maptab[id2]=exchanger;}
}

// Canalized BF: initialization and remove.
void boolMap_Cana::Cana_Initial(){
    boolMap::Initial();
    Selected=(int *)malloc(k*sizeof(int));
    Canalzing=(int *)malloc(k*sizeof(int));
    Canalized=(int *)malloc(k*sizeof(int));
}
void boolMap_Cana::Cana_Reset(){
    boolMap::Reset();
    free(Selected);free(Canalized);free(Canalzing);
    Selected=Canalzing=Canalized=NULL;
}
// Canalized Pure BF: generator.
void boolMap_Cana::Cana_Pure_Gen(){
    boolMap::RBF_p(0.5);// Set initial random mapping table.
    int ii,tmp,*ipt;
    if(Fixed==0||Fixed==1){// Should be set as same.
        Canalized[0]=Canalzing[0]=Fixed;}
    else {// The canalized/canalyzing can be random.
        Canalized[0]=r_binary(mt);
        Canalzing[0]=r_binary(mt);}
    tmp=(int)(k*runif(mt));
    Selected[0]=tmp;
    ipt=enum_dict+tmp*length;
    for(ii=0;ii<length;ii++,ipt++){
        if(*ipt==Canalzing[0])
            maptab[ii]=Canalized[0];}
}
// Canalized Netsed BF: generator.
void boolMap_Cana::Cana_Nest_Gen(){
    boolMap::RBF_p(0.5);// Set initial random mapping table.
    int ii,jj,*ipt,Label[k];
    if(Fixed==0||Fixed==1){// Should be set as same.
        for(ii=0;ii<Deep;ii++){
            Canalized[ii]=Canalzing[ii]=Fixed;}}
    else {// The canalized/canalyzing can be random.
        for(ii=0;ii<Deep;ii++){
            Canalized[ii]=r_binary(mt);
            Canalzing[ii]=r_binary(mt);}}
    memset(Label,0,k*sizeof(int));
    Choose_NK(k,Deep,Selected,Label,-999);// -999: any elements can be chosen.
    for(ii=Deep-1;ii>=0;ii--){// Here need enum_dict, can inverse fill in.
        ipt=enum_dict+length*Selected[ii];// Obtain the index address.
        for(jj=0;jj<length;jj++,ipt++){
            if(*ipt==Canalzing[ii]){
                *(maptab+jj)=Canalized[ii];}}}
}
// General Canalized BF: generator.
void boolMap_Cana::Cana_Gen(){
    switch(Deep){
        case 1:boolMap_Cana::Cana_Pure_Gen();break;// Argu1: Fixed, Argu2: random
        case 2:boolMap_Cana::Cana_Nest_Gen();break;// Argu1: 
        default:printf("Error Type input.\n");getchar();break;}
}
// Canalized BF: identify.
int boolMap_Cana::is_Cana(){// [NOTE] We only consider outermost layer.
    int results=0,logi;
    int ii,jj,kk,*ipt1,*ipt2,boundary_l=1,boundary_u=length>>1;
    int cana_0,cana_1,flag_0=1,flag_1=1;
    for(ii=0;ii<k;ii++){
        ipt1=maptab;
        ipt2=maptab+boundary_l;
        cana_0=*ipt1;// Set 0-input canalized values.
        cana_1=*ipt2;// Set 1-input canalized values.
        logi=flag_0=flag_1=0;
        for(jj=0;jj<boundary_u;jj++){
            for(kk=0;kk<boundary_l;kk++){
                // Keep same, either in-0 or in-1.
                flag_0+=(ipt1[kk]!=cana_0);
                flag_1+=(ipt2[kk]!=cana_1);
                if(flag_0&&flag_1){// Both 0/1 input are not canalized.
                    logi=1;break;}}
            if(logi)break;
            ipt1+=(boundary_l<<1);
            ipt2+=(boundary_l<<1);}
        if(logi){
            boundary_l=boundary_l<<1;
            boundary_u=boundary_u>>1;}
        else {
            results=1;break;}}
    return results;
}

// Post2 BF: generator.
void boolMap_Post::Post_Gen(){
    // [NOTE] Only consider the first element as seed, successor must be intersect each other.
    int ii,jj,tmp,flag,intersect,Selected,Slots[length],*ipt;
    do{// Get first element, intersecting seed, not zero.
        tmp=(int)(length*runif(mt));
        }while(tmp==0);
    Selected=1;labels[tmp]=1;Slots[0]=tmp;
    for(ii=0;ii<(length>>1);ii++){// Only get a loop with 2^k/2 length, regardless of whether code is suitable.
        do{
            tmp=(int)(length*runif(mt));
            }while(labels[tmp]);// Has been selected.
        flag=1;intersect=tmp;// Get seed.
        for(jj=0;jj<Selected;jj++){
            intersect&=Slots[jj];
            if(intersect==0){
                flag=0;break;}}
        if(flag){// Find a suitable slot.
            labels[tmp]=1;
            Slots[Selected++]=tmp;}}
    // Need to judge {A^2} or {a^2}.
    if(Post_Type==1||Post_Type==0){;}// {A^2} or {a^2} class.
    else {// Not set Post_Type.
        Post_Type=r_binary(mt);}
    if(Post_Type==1){// {A^2} subclass
        ipt=maptab;
        for(ii=0;ii<length;ii++,ipt++)*ipt=0;
        // Set "one" position.
        for(ii=0;ii<Selected;ii++){
            labels[Slots[ii]]=0;// Recover label slot.
            maptab[Slots[ii]]=1;}}
    else {// {a^2} subclass
        ipt=maptab;tmp=length-1;
        for(ii=0;ii<length;ii++,ipt++)*ipt=1;
        for(ii=0;ii<Selected;ii++){
            labels[Slots[ii]]=0;// Recover label slot.
            maptab[(~Slots[ii])&tmp]=0;}}// Should get inverse bits.
}
// Post2 BF: identify.
int boolMap_Post::is_Post(){
    int ii,results=1,flag_0,flag_1,*ipt=maptab;
    flag_0=flag_1=(1<<k)-1;
    for(ii=0;ii<length;ii++,ipt++){
        if(*ipt==1){// Map to 1, exist x_i=1.
            flag_1&=ii;}
        else {// Map to 0, exist x_i=0.
            flag_0&=(~ii);}
        if(flag_0+flag_1==0){
            results=0;break;}}
    return results;
}

// Dominant BF: prepare major state.
void boolMap_Domin::Domin_Pre_Major(){
    int ii,jj,num[2],*ipt;
    for(ii=0;ii<length;ii++){
        num[0]=num[1]=0;
        ipt=enum_dict+ii;
        for(jj=0;jj<k;jj++,ipt+=length){
            num[*ipt]++;}
        if(num[0]>num[1]){labels[ii]=1;}
        else if(num[0]<num[1]){labels[ii]=2;}
        else {labels[ii]=3;}}
}
// Dominant BF: generator.
void boolMap_Domin::Domin_Gen(){
    int ii=0,*ipt=labels;
    if(Domin_Type==0||Domin_Type==1){// Domin_Type has been set.
        for(ii=0;ii<length;ii++,ipt++){
            switch(*ipt){
                case 1:maptab[ii]=0;break;
                case 2:maptab[ii]=1;break;
                case 3:maptab[ii]=Domin_Type;break;
                default: printf("Error in Block_Gen(), invalid Para.\n");getchar();}}}
    else {// Has been set.
        for(ii=0;ii<length;ii++,ipt++){
            switch(*ipt){
                case 1:maptab[ii]=0;break;
                case 2:maptab[ii]=1;break;
                case 3:maptab[ii]=r_binary(mt);break;
                default: printf("Error in Block_Gen(), invalid Para.\n");getchar();}}}
}
// Dominant BF: identify.
int boolMap_Domin::is_Domin(){
    int ii,logi=1;
    int *ipt1=maptab,*ipt2=labels;
    for(ii=0;ii<length;ii++){
        if(((*ipt1+1)&(*ipt2))==0){// The state not meet the state-major condition.
            logi=0;break;}
        ipt1++;ipt2++;}
    return logi;
}

// Effective BF: identify.
int boolMap_Effec::is_Effec(){
    int results=0,logi=0;
    int ii,jj,kk,*ipt1,*ipt2,boundary_l=1,boundary_u=length>>1;
    for(ii=0;ii<k;ii++){
        ipt1=maptab;
        ipt2=maptab+boundary_l;
        logi=0;
        for(jj=0;jj<boundary_u;jj++){
            for(kk=0;kk<boundary_l;kk++){
                if((ipt1[kk])!=(ipt2[kk])){// Just one case can meet condition.
                    logi=1;break;}}
            if(logi)break;
            ipt1+=(boundary_l<<1);
            ipt2+=(boundary_l<<1);}
        if(logi){// Meet condition
            boundary_l=boundary_l<<1;
            boundary_u=boundary_u>>1;
            results++;}
        else break;}
    return (results==k);
}
// Effective BF: generator.
void boolMap_Effec::Effec_Gen(){
    int logi=0;
    while(logi==0){
        boolMap::RBF_p(0.50);// Just use random generator.
        logi=boolMap_Effec::is_Effec();}
}

// Monotone BF: identify.
int boolMap_Unate::is_Unate(){
    int ii,jj,kk,logi0=1,logi1=1;
    int exit=0,lower=1,upper=length>>1;
    int *ipt0,*ipt1;
    for(ii=0;ii<k;ii++){// Should ensure all inputs are SAME monotone type.
        ipt0=maptab;
        ipt1=maptab+lower;
        for(jj=0;jj<upper;jj++){
            for(kk=0;kk<lower;kk++){
                if(ipt0[kk]>ipt1[kk])logi1=0;// Contradiction to increase
                if(ipt0[kk]<ipt1[kk])logi0=0;// Contradiction to decrease
                if(logi1+logi0==0){
                    exit=1;// Exit repeat judgment, not increase/decrease monotone.
                    break;}}
            if(exit)break;
            ipt0+=(lower<<1);
            ipt1+=(lower<<1);}
        if(exit)break;
        lower<<=1;// Shrinking range of analysis.
        upper>>=1;}
    return (logi0+logi1>0);
}
int boolMap_Unate::is_Mixed_Unate(){
    int ii,jj,kk,logi0=1,logi1=1;
    int exit=0,lower=1,upper=length>>1;
    int *ipt0,*ipt1;
    for(ii=0;ii<k;ii++){// Each variable can exhibit different monotone type.
        ipt0=maptab;
        ipt1=maptab+lower;
        logi0=logi1=1;
        for(jj=0;jj<upper;jj++){
            for(kk=0;kk<lower;kk++){
                if(ipt0[kk]>ipt1[kk])logi1=0;
                if(ipt0[kk]<ipt1[kk])logi0=0;
                if(logi1+logi0==0){
                    exit=1;
                    break;}}
            if(exit)break;
            ipt0+=(lower<<1);
            ipt1+=(lower<<1);}
        if(exit)break;
        lower<<=1;
        upper>>=1;}
    return (logi0+logi1>0);
}
// Monotone BF: generator.
void boolMap_Unate::Unate_Gen(){
    /*  Set some random seed vectors, and set their successor vectors as candidate.
        Generate an increasing patterns, if (0110010) as one candidate map to one, then all successors should be configurated as one, like (1110010)=>1, (0111010)=>1, ..., (1111111)=>1 (totally 2^4-1=15 successors). */
    int ii,jj,tmp,inverse,n_seed[k]; 
    for(ii=0;ii<k;ii++){
        n_seed[ii]=(int)((length-1)*runif(mt))+1;}// Candidate: 1 ~ length-1, can be repeated.
    memset(maptab,0,length*sizeof(int));
    for(ii=0;ii<k;ii++){
        tmp=n_seed[ii];
        maptab[tmp]=1;
        for(jj=tmp+1;jj<length;jj++){// Forward analysis.
            if((jj&tmp)==tmp){
                maptab[jj]=1;}}}
    if(Unate_Type==0||Unate_Type==1){// Have been set.
        inverse=Unate_Type;}
    else inverse=r_binary(mt);// Have to randomly set one.
    if(inverse==0){// Inverse the relation.
        for(ii=0;ii<length;ii++){
            maptab[ii]=!maptab[ii];}}
}

// Threshold BF: identify.[WARNING] require Z3 solver library
int boolMap_Thres::is_Thres(){
    int ii,jj,res=0;
    z3::context TT_V,TT_W;// Set z3's context.
    std::vector <z3::expr> vars,wars;
    for(ii=0;ii<k;ii++){// Set variable number of VAR slots.
        vars.push_back(TT_V.int_const(("a_" + std::to_string(ii)).c_str()));
        wars.push_back(TT_W.int_const(("b_" + std::to_string(ii)).c_str()));}
    // Build solver.
    z3::expr theta=TT_V.int_const("theta");
    z3::expr exprss1=theta;
    z3::solver solve_it1(TT_V);
    z3::expr gamma=TT_W.int_const("gamma");
    z3::expr exprss2=gamma;
    z3::solver solve_it2(TT_W);
    for(ii=0;ii<length;ii++){
        exprss1=theta;exprss2=gamma;
        for(jj=0;jj<k;jj++){
            if((ii&(1<<jj))){// This bit is one (effective).
                exprss1=exprss1+vars[jj];
                exprss2=exprss2+wars[jj];}}
        if(maptab[ii]>0){
            solve_it1.add(exprss1>=0);
            solve_it2.add(exprss2>0);}
        else {
            solve_it1.add(exprss1<0);
            solve_it2.add(exprss2<=0);}
        for(jj=0;jj<k;jj++){// Regulatory weight value are limited to ±1 (a_{xy}=±1).
            solve_it1.add(1==vars[jj]||-1==vars[jj]);
            solve_it2.add(1==wars[jj]||-1==wars[jj]);}
        solve_it1.add(0==theta);
        solve_it2.add(0==gamma);}// Thershold is zero.
    if(solve_it1.check()==z3::sat){// Zero-summation is set as 1.
        z3::model Models=solve_it1.get_model();
        //for(ii=0;ii<k;ii++){printf("%ld,",Models.eval(vars[ii]).get_numeral_int64());}
        res=1;}
    else if(solve_it2.check()==z3::sat){// Zero-summation is set as 0.
        z3::model Models=solve_it2.get_model();
        //for(ii=0;ii<k;ii++){printf("%ld,",Models.eval(wars[ii]).get_numeral_int64());}
        res=1;}
    else ;
    return res;
}
// Threshold BF: generator.
void boolMap_Thres::Thres_Gen(double p_bias){// [NOTE] Bias p can regulate positive/nagetive inputs.
    int ii,jj,sum=0,rand_fill,pn_in[k],*ipt;
    for(ii=0;ii<k;ii++){// Set +/- inputs. Here is the core of threshold Boolean function.
        pn_in[ii]=((runif(mt)<p_bias)<<1)-1;
        sum+=pn_in[ii];}// Balance the positive/nagetive inputs.
    if(Thres_Type==0||Thres_Type==1){
        rand_fill=Thres_Type;}
    else {
        if(sum!=0){
            rand_fill=(sum<0);}
        else rand_fill=r_binary(mt);}
    for(ii=0;ii<length;ii++){
        sum=0;
        ipt=enum_dict+ii;
        for(jj=0;jj<k;jj++){
            sum+=pn_in[jj]*(*ipt);
            ipt+=length;}
        if(sum>0){
            maptab[ii]=1;}
        else if(sum<0){
            maptab[ii]=0;}
        else maptab[ii]=rand_fill;}
}
#endif // BOOLEAN_FUNCTION_H_INCLUDED
