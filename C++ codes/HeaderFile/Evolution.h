#ifndef EVOLUTION_H_INCLUDED
#define EVOLUTION_H_INCLUDED
#include "Analysis_BioNets_Linux.h"
// Log20230712: merge class BioNets_Random and Evolution to obtain a more easy analyzing code.
// Analysis of real biological networks.
int HammingDistance(int *v1,int *v2,int total){
    int sum=0;
    for(int ii=0;ii<total;ii++){
        sum+=(v1[ii]!=v2[ii]);}
    return sum;
}
int IsSame2Vector(int *v1,int *v2,int total_va,int *valid){
    int logi=1;
    for(int ii=0;ii<total_va;ii++){
        if(v1[valid[ii]]!=v2[valid[ii]]){logi=0;break;}}
    return logi;
}
int UpdateOneNode(int *ss,int deg,int *pars,int *map){
    int ii,sum=0;
    for(ii=0;ii<deg;ii++){
        sum+=(ss[pars[ii]]<<ii);}
    return map[sum];
}
void Dec2Binary(int *ss,int Total,int DecNum){
    int ii,code=DecNum;
    for(ii=0;ii<Total;ii++){
        ss[ii]=code&1;
        code=code>>1;}
}
void SwapSmallestState(int **Looper,int num,int length){// Let smallest be the first, to match parallel results.
    int ii,jj,logi,min_id=0,tmp[length],*ipt,bits=length*sizeof(int);
    for(ii=1;ii<num;ii++){// memcpy(tmp,Looper[0],bits);
        ipt=Looper[min_id];
        for(jj=0;jj<length;jj++){
            if(ipt[jj]>Looper[ii][jj]){
                min_id=ii;
                break;}}}
    if(min_id!=0){// Changed!
        ipt=Looper[0];Looper[0]=Looper[min_id];Looper[min_id]=ipt;}
}
void FreeStatePrint(FILE *fp,int *ss,int total,int *id_slot){
    fprintf(fp,"%1d",ss[id_slot[0]]);
    for(int ii=1;ii<total;ii++)fprintf(fp,"\t%1d",ss[id_slot[ii]]);
    fprintf(fp,"\n");
}
// Sub_Class: random assignment 
class BioNets_Random:public BioNets{
    public:
        double Av_Deg;
        void Generator(int total,int init_total,double av_deg);
        int *StateVector;
        int *Candidate;
        int ***Loops;
        int *loop_num;
        int valid_num;// Not all nodes are valid.
        int *valid_code;// Record valid node's code.
        void PrintGene(FILE *fp);
        void PrintState(FILE *fp);
        void FindAttractorLimitcycle(FILE *fp1,FILE *fp2,int TimeNum,int Validss);
        void sub1_ValidNodeCode();
    private:
        int Update();
        void RandomVec();
        int Attractor();
        int LimitCyclePool_NewLoop(int index);
        int LimitCyclePool(int *RecordVec,int MaxLeng,int *GroupSize);
        int sub2_InverseCheck(int WholeID);
};
void BioNets_Random::Generator(int total,int init_total,double av_deg){
    Av_Deg=av_deg;
    Gene_Network.total=total; 
    GeneNames=(char**)malloc(total*sizeof(char*));
    Maps=(int**)malloc(total*sizeof(int*));
    Parents=(int **)malloc(total*sizeof(int *));
    Gene_Network.Initial_Network();
    int ii,jj,kk,tmp,*ind=Gene_Network.InDeg,*otd=Gene_Network.OtDeg;
    struct Regulation *net=Gene_Network.Network;
    // Build a directed ER network.
    std::poisson_distribution<int> PoissonDist((int)Av_Deg);
    // Assign candiate momery.
    std::vector<int> id_num;//id_num.resize(init_total);
    //for(ii=0;ii<init_total;ii++)id_num[ii]=ii;
    for(ii=0;ii<init_total;ii++)id_num.push_back(ii);
    // Part1: assigned linking nodes.
    for(ii=0;ii<init_total;ii++){
        tmp=PoissonDist(mt);
        if(0==tmp)tmp=1;else if(tmp>=init_total)tmp=init_total-1;else ;
        GeneNames[ii]=(char*)malloc(16*sizeof(char));
        Parents[ii]=(int *)malloc(tmp*sizeof(int));
        Maps[ii]=(int *)malloc((1<<tmp)*sizeof(int));
        snprintf(GeneNames[ii],16,"G%03d",ii);
        ind[ii]=tmp;
        std::shuffle(id_num.begin(),id_num.end(),mt);// Shuffle the vector.
        jj=0;kk=0;
        while(kk<tmp){
            if(id_num[jj]!=ii){
                Parents[ii][kk]=id_num[jj];
                otd[id_num[jj]]++;
                BuildRegulationship(net,ii,id_num[jj]);
                kk++;}
            jj++;}
        for(jj=0;jj<(1<<tmp);jj++){
            Maps[ii][jj]=runif(mt)>0.5;}}
    // Part2: null nodes of system (for later evoluting).
    for(ii=init_total;ii<total;ii++){
        GeneNames[ii]=(char*)malloc(16*sizeof(char));
        //sprintf_s(GeneNames[ii],16,"G%03d",ii);
        snprintf(GeneNames[ii],16,"G%03d",ii);
        Parents[ii]=Maps[ii]=NULL;}
}
// Auxiliary functions of attractor pools.
void BioNets_Random::RandomVec(){
    //for(int ii=0;ii<Gene_Network.total;ii++)StateVector[ii]=runif(mt)>0.5;
    for(int ii=0;ii<valid_num;ii++)StateVector[valid_code[ii]]=runif(mt)>0.5;
}
void BioNets_Random::PrintGene(FILE *fp){
    fprintf(fp,"%s",GeneNames[valid_code[0]]);
    for(int ii=1;ii<valid_num;ii++)fprintf(fp,"\t%s",GeneNames[valid_code[ii]]);
    fprintf(fp,"\n");
}
void BioNets_Random::PrintState(FILE *fp){
    fprintf(fp,"%1d",StateVector[valid_code[0]]);
    for(int ii=1;ii<valid_num;ii++)fprintf(fp,"\t%1d",StateVector[valid_code[ii]]);
    fprintf(fp,"\n");
}
int BioNets_Random::Update(){
    int total=valid_num,*degs=Gene_Network.InDeg;
    int ii,jj=0,rand_id,tmps,tmp_int[total];
    for(ii=0;ii<valid_num;ii++){
        if(Candidate[valid_code[ii]]){
            tmps=UpdateOneNode(StateVector,degs[valid_code[ii]],Parents[valid_code[ii]],Maps[valid_code[ii]]);
            if(tmps==StateVector[valid_code[ii]]){
                Candidate[valid_code[ii]]=0;}
            else {
            tmp_int[jj++]=valid_code[ii];}}}
    if(jj>0){
        rand_id=(int)(runif(mt)*jj);
        int tmps;
        rand_id=tmp_int[rand_id];// Random select one, 
        StateVector[rand_id]=UpdateOneNode(StateVector,degs[rand_id],Parents[rand_id],Maps[rand_id]);
        struct Regulation *pt=Gene_Network.Network+rand_id;
        Candidate[rand_id]=0;
        while(pt->next!=NULL){
            tmps=sub2_InverseCheck(pt->next->code);
            Candidate[valid_code[tmps]]=1;
            pt=pt->next;}}
    return jj;// "jj==0" means has reached an attractor.
}
int BioNets_Random::Attractor(){
    int changes=1,step=0,MaxStep=1000;// MAX step:1000
    while(changes>0&&step<MaxStep){
        changes=BioNets_Random::Update();
        step++;}
    return (step<MaxStep);// "Step<Max", succeed finding.
}
int BioNets_Random::LimitCyclePool_NewLoop(int index){
    int total=Gene_Network.total,*degs=Gene_Network.InDeg;
    int number=1,**iptt,**NewLoop;
    int ii,jj=0,rand_id,tmps,logi_1,tmp_int[total],MaxPool=1024; int intermediate;// for partial network.
    struct Regulation *rpt2,*rpt,*source=(struct Regulation *)malloc(sizeof(struct Regulation));
    source->code=0;source->prev=source->next=NULL;
    NewLoop=(int**)malloc(sizeof(int*));
    NewLoop[0]=(int*)malloc(total*sizeof(int));
    memcpy(NewLoop[0],StateVector,total*sizeof(int));
    for(ii=0;ii<valid_num;ii++){
        intermediate=valid_code[ii];
        if(degs[intermediate]){
            tmps=UpdateOneNode(StateVector,degs[intermediate],Parents[intermediate],Maps[intermediate]);
            if(tmps!=StateVector[intermediate]){// Add candiate SS.
                iptt=NewLoop;
                NewLoop=(int**)malloc((number+1)*sizeof(int*));
                memcpy(NewLoop,iptt,number*sizeof(int*));free(iptt);
                NewLoop[number]=(int*)malloc(total*sizeof(int));// Record new SS.
                memcpy(NewLoop[number],StateVector,total*sizeof(int));
                NewLoop[number][ii]=tmps;
                rpt=(struct Regulation *)malloc(sizeof(struct Regulation));// Save tree structure.
                rpt->next=NULL;//ppp++;
                rpt->code=number;
                if(source->next==NULL){
                    source->next=rpt;}
                else {
                    rpt->next=source->next;
                    source->next=rpt;}
                number++;}}
    }
    while(source->next!=NULL&&number<MaxPool){
        memcpy(StateVector,NewLoop[source->next->code],total*sizeof(int));
        rpt=source->next;source->next=source->next->next;free(rpt);//nnn++;
        memcpy(tmp_int,StateVector,total*sizeof(int));
        for(ii=0;ii<valid_num;ii++){
            intermediate=valid_code[ii];
            if(degs[intermediate]&&number<MaxPool){
                tmps=UpdateOneNode(StateVector,degs[intermediate],Parents[intermediate],Maps[intermediate]);
                if(tmps!=StateVector[intermediate]){// Add candiate SS.
                    tmp_int[intermediate]=tmps;
                    logi_1=1;
                    for(jj=0;jj<number;jj++){
                        if(IsSame2Vector(NewLoop[jj],tmp_int,valid_num,valid_code)){logi_1=0;break;}}
                    if(logi_1){// Not existed.
                        iptt=NewLoop;
                        NewLoop=(int**)malloc((number+1)*sizeof(int*));
                        memcpy(NewLoop,iptt,number*sizeof(int*));free(iptt);
                        NewLoop[number]=(int*)malloc(total*sizeof(int));// Record new SS.
                        memcpy(NewLoop[number],StateVector,total*sizeof(int));
                        NewLoop[number][ii]=tmps;
                        rpt=(struct Regulation *)malloc(sizeof(struct Regulation));// Save tree structure.
                        rpt->next=NULL;
                        rpt->code=number;//ppp++;
                        if(source->next==NULL){
                            source->next=rpt;}
                        else {
                            rpt->next=source->next;
                            source->next=rpt;}
                        number++;}}
                tmp_int[intermediate]=StateVector[intermediate];}}}
    if(number<MaxPool){// Small loop <1024
        free(source);
        SwapSmallestState(NewLoop,number,total);// Let smallest is the head.
        *(Loops+index)=NewLoop;}
    else {// Too large, should be carefully.
        *(Loops+index)=NULL;
        for(ii=0;ii<number;ii++){
            free(NewLoop[ii]);}
        while(source->next!=NULL){
            rpt=source->next;source->next=source->next->next;free(rpt);}
        free(NewLoop);free(source);
        number=-number;}// Minus denotes a large loop.
    return number;
}
int BioNets_Random::LimitCyclePool(int *RecordVec,int MaxLeng,int *GroupSize){
    int ii,jj,dist,changes,logi,linshi;
    int ***ipttt,*ipt;int which_group;
    //FILE *fpt=fopen("Y://POOL.txt","w");
    if(loop_num[0]==0){// The first loop-pool, small or large should be recorded.
        Loops=(int***)malloc(sizeof(int**));
        ipt=loop_num;
        loop_num=(int*)malloc(2*sizeof(int));
        memcpy(loop_num,ipt,sizeof(int));
        loop_num[1]=LimitCyclePool_NewLoop(loop_num[0]);
        loop_num[0]++;
        which_group=0;
        *GroupSize=loop_num[which_group+1];}
    else {
        logi=1;
        for(ii=0;ii<loop_num[0];ii++){
            if(logi){
                for(jj=0;jj<loop_num[ii+1];jj++){
                    if(IsSame2Vector(Loops[ii][jj],RecordVec,valid_num,valid_code)){
                        which_group=ii;*GroupSize=loop_num[ii+1];
                        logi=0;break;}}}}
        if(logi==1){// New complex loop.
            ipttt=Loops;//printf("Has a new pool.\n");
            Loops=(int***)malloc((loop_num[0]+1)*sizeof(int**));
            memcpy(Loops,ipttt,sizeof(int**)*loop_num[0]);free(ipttt);
            linshi=BioNets_Random::LimitCyclePool_NewLoop(loop_num[0]);
            if(-1024==linshi){// Reject! Maybe a super-large pool.
                // Return a small slot;
                for(ii=0;ii<loop_num[0];ii++){
                    if(-1024==loop_num[ii+1]){
                        break;}}
                if(ii>=loop_num[0]){// Not existed!!
                    ipt=loop_num;
                    loop_num=(int*)malloc((ipt[0]+2)*sizeof(int));
                    memcpy(loop_num,ipt,(ipt[0]+1)*(sizeof(int)));free(ipt);
                    loop_num[loop_num[0]+1]=linshi;
                    which_group=loop_num[0];
                    loop_num[0]++;}
                else {//Has existed!!
                    ipttt=Loops;//printf("Has a new pool.\n");
                    Loops=(int***)malloc((loop_num[0])*sizeof(int**));
                    memcpy(Loops,ipttt,sizeof(int**)*loop_num[0]);free(ipttt);
                    which_group=ii;}
                *GroupSize=-1024;}
            else {// Accept! A small pool.
                ipt=loop_num;// printf("Accept one.\n");
                loop_num=(int*)malloc((ipt[0]+2)*sizeof(int));
                memcpy(loop_num,ipt,(ipt[0]+1)*(sizeof(int)));free(ipt);
                loop_num[loop_num[0]+1]=linshi;
                which_group=loop_num[0];
                loop_num[0]++;
                *GroupSize=loop_num[loop_num[0]];}}
        else ;}
    return which_group;// Regradless of repeated, record them.
}
void BioNets_Random::FindAttractorLimitcycle(FILE *fp1,FILE *fp2,int TimeNum,int Validss){
    int total=Gene_Network.total,Flag[total],Times;
    int loop_type,size[]={1};
    StateVector=(int*)malloc(total*sizeof(int));
    Candidate=(int*)malloc(total*sizeof(int));
    Loops=NULL;loop_num=(int *)malloc(sizeof(int));loop_num[0]=0;// Save each loop number.
    //BioNets_Random::PrintGene(fp1);
    fprintf(fp2,"period\n");
    if(valid_num<24){// If system size is less than 24.
        Times=1<<valid_num;
        for(int ii=0;ii<Times;ii++){
            memcpy(Candidate,Gene_Network.InDeg,total*sizeof(int));
            BioNets_Random::RandomVec();
            if(BioNets_Random::Attractor()){// Is attractor, no need to saved.
                fprintf(fp2,"1\n");}
            else {
                memcpy(Flag,StateVector,total*sizeof(int));
                loop_type=BioNets_Random::LimitCyclePool(Flag,1000,size);
                fprintf(fp2,"%d\n",*size);}}}
    else {
        Times=TimeNum;
        for(int ii=0;ii<Times;ii++){
            memcpy(Candidate,Gene_Network.InDeg,total*sizeof(int));
            BioNets_Random::RandomVec();
            if(BioNets_Random::Attractor()){// Is attractor.
                fprintf(fp2,"1\n");}
            else {
                memcpy(Flag,StateVector,total*sizeof(int));
                loop_type=BioNets_Random::LimitCyclePool(Flag,1000,size);
                fprintf(fp2,"%d\n",*size);}}}
    for(int ii=0;ii<loop_num[0];ii++){
        for(int jj=0;jj<loop_num[ii+1];jj++){
            free(Loops[ii][jj]);}
        free(Loops[ii]);}
    free(StateVector);free(Candidate);free(Loops);
    free(loop_num); 
}

void BioNets_Random::sub1_ValidNodeCode(){
    int ii,sum=0,total=Gene_Network.total,*ind=Gene_Network.InDeg,*otd=Gene_Network.OtDeg;
    free(valid_code);
    valid_code=(int*)malloc(valid_num*sizeof(int));
    for(ii=0;ii<total;ii++){
        if(ind[ii]||otd[ii]){// Node is not isolated!
            valid_code[sum]=ii;
            sum++;}}
    if(sum!=valid_num){
        printf("Error in obtaining valid coder!\n");getchar();}
}
int BioNets_Random::sub2_InverseCheck(int WholeID){
    int ids=-666;//[NOTE] the whole-ID must be included in valid_code slot.
    for(int ii=0;ii<valid_num;ii++){
        if(WholeID==valid_code[ii]){
            ids=ii;break;}}
    return ids;
}

// Independent class: evolution.
class Evolution:public BioNets_Random{
    public:
        //int UpperNode;
        int ValidNode;
        // Evolution of a biological network.
        int Evo_Patterns_Rand(double *Prob);
        
    private:
        void Evo_DelNode_Maps_Pars(int target,int source,int *indeg);
        void Evo_AddNode(int Type);
        void Evo_DelNode();
        void Evo_SwapEdge();
        void Evo_ChangeMap(int Type,int Inverse);
        void Evo_AddEdge(int NewNode_Type);
        void Evo_DelEdge();
        int IsolatedNode();
        int RandomChooseWays(double *prob);
};
int Evolution::IsolatedNode(){
    int sum=0,*ind=Gene_Network.InDeg,*otd=Gene_Network.OtDeg;
    for(int ii=0;ii<Gene_Network.total;ii++){
        if(ind[ii]||otd[ii])sum++;}
    return sum;
}
// Add a new node into system.
void Evolution::Evo_AddNode(int Type){// 0-random,1-cana, 2-thre, minus-random
    int ID=-666,tmp,tmp2,type=Type,*indeg=Gene_Network.InDeg,*otdeg=Gene_Network.OtDeg;
    for(tmp=0;tmp<Gene_Network.total;tmp++){
        if(indeg[tmp]==0&&otdeg[tmp]==0){
            ID=tmp;break;}}
    // How many input?
    if(ID==-666){
        printf("Empty System!\n");}
    std::poisson_distribution<int> PoissonDist((int)Av_Deg);
    tmp=PoissonDist(mt);
    if(0==tmp)tmp=1;if(tmp>=ValidNode)tmp=ValidNode-1;
    if(tmp<0){printf("Abnormal inputs");}
    Parents[ID]=(int*)malloc(tmp*sizeof(int));
    Maps[ID]=(int*)malloc((1<<tmp)*sizeof(int));
    // Set parents
    std::vector<int> candiate;
    for(tmp2=0;tmp2<Gene_Network.total;tmp2++){
        if(indeg[tmp2]!=0||otdeg[tmp2]!=0){// Record all non-in/ot-degree nodes.
            candiate.push_back(tmp2);}}
    indeg[ID]=tmp;
    std::random_shuffle(candiate.begin(),candiate.end());
    for(tmp=0;tmp<indeg[ID];tmp++){
        Parents[ID][tmp]=candiate[tmp];
        otdeg[candiate[tmp]]++;
        BuildRegulationship(Gene_Network.Network,ID,candiate[tmp]);}
    // Set a mapping table [Only Cana/Thres]
    boolMap_Cana bf_c;
    boolMap_Thres bf_t;
    boolMap bf_r;
    if(Type<0){
        type=(int)(3.0*runif(mt));}
    if(type==1){
        bf_c.k=indeg[ID]; bf_c.length=1<<indeg[ID];
        bf_c.Deep=1; bf_c.Fixed=9;// 9: random canalized and canalizing.
        bf_c.Cana_Initial(); bf_c.In_Enumerate();
        bf_c.maptab=Maps[ID];
        bf_c.Cana_Gen();
        bf_c.Cana_Reset();}
    else if(type==2){
        bf_t.k=indeg[ID]; bf_t.length=1<<indeg[ID];
        bf_t.Initial(); bf_t.In_Enumerate();
        bf_t.maptab=Maps[ID];
        bf_t.Thres_Gen(0.500);
        bf_t.Reset();}
    else if(type==0){// totally randomly 
        bf_r.k=indeg[ID]; bf_r.length=1<<indeg[ID];
        bf_r.maptab=Maps[ID];
        bf_r.RBF_p(0.5);}//bf_r.Reset();
    else {
        printf("Error input.\n");getchar();}
}
// Randomly delete a node.
void Evolution::Evo_DelNode_Maps_Pars(int target,int source,int *indeg){
    int loc,*ipt,tmp1=0,tmp2=0;//[NOTE] Some indeg=1 nodes issue.
    int ii,jj;
    // Update pars.
    if(indeg[target]>1){
        ipt=Parents[target];
        Parents[target]=(int*)malloc((indeg[target]-1)*sizeof(int));
        for(ii=0;ii<indeg[target];ii++){
            if(ipt[ii]==source){
                loc=ii;tmp1++;}
            else {
                Parents[target][tmp2++]=ipt[tmp1++];}}
        free(ipt);}
    else {
        free(Parents[target]);Parents[target]=NULL;}
    // Update mapping
    if(indeg[target]>1){
        ipt=Maps[target];jj=0;
        tmp1=1<<(indeg[target]-1);
        tmp2=1<<loc;
        Maps[target]=(int*)malloc(tmp1*sizeof(int*));
        for(ii=0;ii<(tmp1<<1);ii++){
            if((ii&tmp2)==0){// bit at loc is zero.
                Maps[target][jj++]=ipt[ii];}}
        free(ipt);}
    else {
        free(Maps[target]);Maps[target]=NULL;}
}
void Evolution::Evo_DelNode(){
    int ID,tmp,*indeg=Gene_Network.InDeg,*otdeg=Gene_Network.OtDeg;
    std::vector<int> candiate;int flag=0;
    //candiate.resize(ValidNode);
    for(tmp=0;tmp<Gene_Network.total;tmp++){
        if(indeg[tmp]!=0||otdeg[tmp]!=0){
            flag++;
            candiate.push_back(tmp);}}
    if(flag!=ValidNode){printf("error-1!!!");getchar();}
    if(0==ValidNode){printf("error-2!!!");getchar();}
    ID=(int)(runif(mt)*ValidNode);
    ID=candiate[ID];
    struct Regulation *rpt,*del;
    rpt=&(Gene_Network.Network[ID]);
    while(rpt->next!=NULL){// Remove child links.
        tmp=rpt->next->code;
        Evo_DelNode_Maps_Pars(tmp,ID,indeg);
        DeleteParents(&Gene_Network.Network[tmp],ID);
        indeg[tmp]--;
        rpt=rpt->next;}
    rpt=&(Gene_Network.Network[ID]);
    while(rpt->prev!=NULL){
        tmp=rpt->prev->code;
        DeleteChilds(&Gene_Network.Network[tmp],ID);
        otdeg[tmp]--;
        rpt=rpt->prev;}
    free(Parents[ID]);free(Maps[ID]);Parents[ID]=Maps[ID]=NULL;// Free pars & maps slot.
    otdeg[ID]=indeg[ID]=0;
    rpt=&(Gene_Network.Network[ID]);// Remove the structure of Evo_Net.network.
    while(rpt->next!=NULL){
        del=rpt->next; rpt->next=del->next; free(del);}
    rpt=&(Gene_Network.Network[ID]);
    while(rpt->prev!=NULL){
        del=rpt->prev; rpt->prev=del->prev; free(del);}
    rpt=&(Gene_Network.Network[ID]);
    if(rpt->next!=NULL||rpt->prev!=NULL){printf("Error Pointers!\n");getchar();}
}
// Swap two edges, A->a + B->b  ==>  A->b + B->a
void Evolution::Evo_SwapEdge(){
    int ii,sou1,sou2,tar1,tar2,maxtest=0,maxtest2,logi=1;
    int *indeg=Gene_Network.InDeg,*otdeg=Gene_Network.OtDeg;
    struct Regulation *rpt,*npt=Gene_Network.Network;
    while(maxtest<50){
        maxtest2=0;
        do{// Search the sources.
            sou1=(int)(runif(mt)*Gene_Network.total);
            sou2=(int)(runif(mt)*Gene_Network.total);
            maxtest2++;
            if(maxtest2>60){logi=0;break;}
        }while(sou1==sou2||otdeg[sou1]==0||otdeg[sou2]==0);
        // Random childs
        if(0==logi){maxtest=666;break;}
        tar1=(int)(runif(mt)*otdeg[sou1]);
        tar2=(int)(runif(mt)*otdeg[sou2]);
        rpt=&npt[sou1];
        for(ii=0;ii<tar1;ii++)rpt=rpt->next;
        tar1=rpt->next->code;
        rpt=&npt[sou2];
        for(ii=0;ii<tar2;ii++)rpt=rpt->next;
        tar2=rpt->next->code;
        if(sou1==tar2||sou2==tar1||ExistChild(&npt[sou1],tar2)||ExistChild(&npt[sou2],tar1))maxtest++;// Find failure!
        else {
            DeleteChilds(&npt[sou1],tar1);DeleteChilds(&npt[sou2],tar2);
            DeleteParents(&npt[tar1],sou1);DeleteParents(&npt[tar2],sou2);
            BuildRegulationship(npt,tar2,sou1);// net, child, parent
            BuildRegulationship(npt,tar1,sou2);// net, child, parent
            break;}}
    if(maxtest>=50)logi=0;
    // Update the parents slots.
    if(logi==1){// Succeed!
        for(ii=0;ii<indeg[tar1];ii++){// tar1<-sou1  ==> tar1<-sou2
            if(Parents[tar1][ii]==sou1){
                Parents[tar1][ii]=sou2;break;}}
        for(ii=0;ii<indeg[tar2];ii++){// tar2<-sou2  ==> tar2<-sou1
            if(Parents[tar2][ii]==sou2){
                Parents[tar2][ii]=sou1;break;}}}
}
// Changed one node's mapping type by single input variable.
void ReturnPointedBit(int bits,int loc,int *map0,int *map1){
    int ii,num0=0,num1=0;
    int length=1<<bits,logi=1<<loc;
    for(ii=0;ii<length;ii++){
        if(ii&logi)map1[num1++]=ii;
        else map0[num0++]=ii;}
}
void ChangedMap_Cana(int bits,int loc,int *map,int inverse){// "inverse" only ensure the canalizing variable 
    int tmp[1<<bits],reps=1<<(bits-1);
    int *id0=tmp,*id1=tmp+reps,*ipt;
    ReturnPointedBit(bits,loc,id0,id1);
    int cana_in,cana_ot;
    cana_in=runif(mt)>0.500;//printf("<%d>",cana_in);
    if(cana_in)ipt=id1;// Regulting "one"
    else ipt=id0;// Regulating "zero"
    if(inverse){// Keep an opposite value.
        cana_ot=1-cana_in;}
    else {
        cana_ot=cana_in;}// Keep same value.
    for(int ii=0;ii<reps;ii++){
        map[ipt[ii]]=cana_ot;}
}
void ChangedMap_Mono(int bits,int loc,int *map,int inverse){// Only set one variable's mapping, not all variables'.
    int tmp[1<<bits],reps=1<<(bits-1),swaper;
    int *id0=tmp,*id1=tmp+reps;
    ReturnPointedBit(bits,loc,id0,id1);
    // [NOTE] if not meet cerita, only exchange two mapping points.
    if(inverse){// Decreasing
        for(int ii=0;ii<reps;ii++){
            if(map[id0[ii]]<map[id1[ii]]){// 0-part < 1-part
                swaper=map[id0[ii]];
                map[id0[ii]]=map[id1[ii]];
                map[id1[ii]]=swaper;}}}
    else {// Increasing
        for(int ii=0;ii<reps;ii++){
            if(map[id0[ii]]>map[id1[ii]]){// 0-part > 1-part
                swaper=map[id0[ii]];
                map[id0[ii]]=map[id1[ii]];
                map[id1[ii]]=swaper;}}}
}
void Evolution::Evo_ChangeMap(int Type,int Inverse){// type: canalized-1, mono-2, minus-random
    int type,ID,Loc,*indeg=Gene_Network.InDeg;
    do{
        ID=(int)(runif(mt)*Gene_Network.total);}while(indeg[ID]==0);
    Loc=(int)(runif(mt)*indeg[ID]);
    if(Type<0){
        type=(runif(mt)>0.5);}
    else {type=Type;}
    switch(type){
        case 0:ChangedMap_Cana(indeg[ID],Loc,Maps[ID],Inverse);break;
        case 1:ChangedMap_Mono(indeg[ID],Loc,Maps[ID],Inverse);break;
        default:printf("Error Changing Map Type.\n");break;}
}
// Add a new edge between existed valid nodes (Not introduce new nodes).
void Evolution::Evo_AddEdge(int NewNode_Type){
    int tar_id,sou_id,*otd=Gene_Network.OtDeg,*ind=Gene_Network.InDeg;
    struct Regulation *rpt,*npt=Gene_Network.Network;
    int logi=1,counter=0;
    do{
        sou_id=(int)(runif(mt)*Gene_Network.total);
        tar_id=(int)(runif(mt)*Gene_Network.total);
        if(counter>20){logi=0;break;}
        counter++;
    }while(sou_id==tar_id|| ExistChild(&npt[sou_id],tar_id) ||
        (otd[sou_id]==0&&ind[sou_id]==0) || (otd[tar_id]==0&&ind[tar_id]==0) );
    if(logi){
        //printf(" %d --> %d ",sou_id,tar_id);
        BuildRegulationship(npt,tar_id,sou_id);
        int *ipt;
        if(0==ind[tar_id]){// Empty of parnets input.
            Parents[tar_id]=(int*)malloc(sizeof(int));
            Parents[tar_id][0]=sou_id;
            Maps[tar_id]=(int*)malloc(2*sizeof(int));
            Maps[tar_id][0]=runif(mt)>0.5;Maps[tar_id][1]=runif(mt)>0.5;}
        else {
            // Update parents.
            ipt=Parents[tar_id];
            Parents[tar_id]=(int*)malloc((1+ind[tar_id])*sizeof(int));
            memcpy(Parents[tar_id],ipt,ind[tar_id]*sizeof(int));
            Parents[tar_id][ind[tar_id]]=sou_id;free(ipt);
            // Update mapping table.
            int length=1<<ind[tar_id];
            ipt=Maps[tar_id];
            Maps[tar_id]=(int*)malloc(2*length*sizeof(int));
            memcpy(Maps[tar_id],ipt,length*sizeof(int));
            int type_fun=NewNode_Type;
            if(type_fun<0){type_fun=(int)(runif(mt)*3.0);}
            int *ip0=Maps[tar_id],*ip1;ip1=ip0+length;
            if(0==type_fun){// random
                for(int kk=0;kk<length;kk++){
                    ip1[kk]=runif(mt)>0.5;}}
            else if(1==type_fun){// cana-1
                int cana=runif(mt)>0.5;
                for(int kk=0;kk<length;kk++){
                    ip1[kk]=cana;}}
            else if(2==type_fun){// mono-2
                int in_de=runif(mt)>0.5;
                if(in_de){// Increasing
                    for(int kk=0;kk<length;kk++){
                        ip1[kk]=(ip0[kk]+(runif(mt)>0.5))>0;}}
                else {// Decreasing
                    for(int kk=0;kk<length;kk++){
                        ip1[kk]=ip0[kk]*(runif(mt)>0.5);}}}
            for(int ok=0;ok<2*length;ok++){
                if(ip0[ok]==0||ip0[ok]==1);
                else {printf("Mapping error.\n");}}}
        otd[sou_id]++;
        ind[tar_id]++;}
}
// Delete an existed edge.
void Evolution::Evo_DelEdge(){
    int ii,ID,tmp,*ipt,*otd=Gene_Network.OtDeg,*ind=Gene_Network.InDeg;
    int target_id,source_id,parent_id,counter=0,logi=1;
    struct Regulation *rpt,*npt=Gene_Network.Network;
    while(1){// Delete an edge randomly.
        source_id=(int)(runif(mt)*Gene_Network.total);counter++;
        if(counter>20){logi=0;break;}
        if(otd[source_id])break;}
    if(logi){
        tmp=(int)(runif(mt)*otd[source_id]);
    rpt=&npt[source_id];
    for(ii=0;ii<tmp;ii++)rpt=rpt->next;
    target_id=rpt->next->code;// Child's code.
    DeleteChilds(&npt[source_id],target_id);
    DeleteParents(&npt[target_id],source_id);
    // Obtain mapping id.
    if(1==ind[target_id]){// Only "source_id" input.
        free(Parents[target_id]);free(Maps[target_id]);
        Parents[target_id]=Maps[target_id]=NULL;}
    else {
        int tmp[1<<ind[target_id]],reps=1<<(ind[target_id]-1),swaper;
        int *id0=tmp,*id1=tmp+reps;
        for(int cor=0;cor<ind[target_id];cor++){
            if(Parents[target_id][cor]==source_id){
                parent_id=cor;break;}}
        ReturnPointedBit(ind[target_id],parent_id,id0,id1);
        ipt=Maps[target_id];
        Maps[target_id]=(int*)malloc(reps*sizeof(int));
        for(int kk=0;kk<reps;kk++){
            Maps[target_id][kk]=ipt[id0[kk]];}
        free(ipt);
        ipt=Parents[target_id];
        Parents[target_id]=(int*)malloc((ind[target_id]-1)*sizeof(int));
        swaper=0;
        for(int cor=0;cor<ind[target_id];cor++){
            if(ipt[cor]!=source_id){
                Parents[target_id][swaper++]=ipt[cor];;}}
        free(ipt);}
    otd[source_id]--;ind[target_id]--;}
    else {;}
}
// Randomly select one type evolution.
int Evolution::RandomChooseWays(double *prob){// Swap, Map, Add, Del
    double ways=runif(mt);
    if(ways<prob[0])return 0;
    else if(ways<prob[1])return 1;
    else if(ways<prob[2])return 2;
    else if(ways<prob[3])return 3;
    else if(ways<prob[4])return 4;
    else return 5;
}
int Evolution::Evo_Patterns_Rand(double *Prob){
    int type=RandomChooseWays(Prob);
    //int type=(int)(6.0*runif(mt));
    if(ValidNode<5&&(type!=2))type=666;
    if((ValidNode>=Gene_Network.total||ValidNode<2)&&type==2)type=666;// Avoid empty system.
    switch(type){
        case 0:Evo_SwapEdge();break;
        case 1:Evo_ChangeMap(-666,0);break;
        case 2:Evo_AddNode((int)(runif(mt)*3)+0);ValidNode++;break;// cana or thres; if consider RANDOM, set it as 0, RANDOM+C/T set it as (int)(runif(mt)*3)
        case 3:Evo_DelNode();ValidNode=IsolatedNode();break;
        case 4:Evo_AddEdge((int)(runif(mt)*2)+1);break;
        case 5:Evo_DelEdge();ValidNode=IsolatedNode();break;
        default:break;}// Jump this loop due to non-meeting conditions.
    return type;
}

void Evoluting_OnlyEdge_Map(char *Dir,char *FileName,int Repeat,double *Probs,int EvoType){
    char target[64],ss[256],tt[512];int linshi;
    FILE *fp_sen_av;
    std::string FilePathName;
    FilePathName.assign(Dir).append("//").append(FileName);
    strcpy(tt,FilePathName.c_str());
    fp_sen_av=fopen(strcat(tt,".txt"),"w");
    Evolution eva;
    for(int ii=0;ii<Repeat;ii++){
        eva.Generator(100,50,5);
        eva.ValidNode=50;
        fprintf(fp_sen_av,"%d,%d,%f,%f\n",eva.ValidNode,eva.TotalEdge(),eva.TotalBias(),eva.AverageSensitivity());
        for(int jj=1;jj<500;jj++){
            //printf("%d,%d,",ii,jj);
            switch(EvoType){
                case 1:linshi=eva.Evo_Patterns_Rand(Probs);break;
                //case 2:linshi=eva.Evo_Patterns_Rand(Probs);break;
                //case 3:linshi=eva.Evo_Patterns_Rand(Probs);break;
                default:printf("Error!\n");getchar();break;}
            //if(eva.Debug()){printf("error-1!");getchar();}
            //if(eva.Debug2()){printf("error-2!");getchar();}
            //if(eva.Debug3()){printf("error-3!");getchar();}
            //printf("%d,%d,%d\n",linshi,eva.ValidNode,eva.TotalEdge());
            fprintf(fp_sen_av,"%d,%d,%f,%f\n",eva.ValidNode,eva.TotalEdge(),eva.TotalBias(),eva.AverageSensitivity());}}
    eva.DeleteBioNetwork();
    fclose(fp_sen_av);
}
void Evoluting_EvoAttrLoop(char *Dir,char *FileName,double *Probs,int Generation){
    char target[64],ss[256],tt[512],num[16];
    FILE *fp_stab=NULL,*fp_step=NULL;
    std::string FilePathName;
    Evolution eva;
    eva.valid_code=NULL;
    eva.Generator(100,50,5);
    eva.ValidNode=eva.valid_num=50;
    // Here evolve the pointed generations.
    for(int jj=1;jj<Generation;jj++){
        eva.Evo_Patterns_Rand(Probs);}
    FilePathName.assign(Dir).append("//").append(FileName).append("_end");//.append(num);
    //strcpy(tt,FilePathName.c_str());
    //fp_stab=fopen(strcat(tt,"_stab.txt"),"w");// No need to record detail attractors.
    strcpy(tt,FilePathName.c_str());
    fp_step=fopen(strcat(tt,"_step.txt"),"w");
    eva.valid_num=eva.ValidNode;
    eva.sub1_ValidNodeCode();// Return valid coders.
    eva.FindAttractorLimitcycle(fp_stab,fp_step,1000000,eva.ValidNode);
    fclose(fp_step);
    free(eva.valid_code);eva.valid_code=NULL;
    eva.DeleteBioNetwork();
}
#endif // EVOLUTION_H_INCLUDED