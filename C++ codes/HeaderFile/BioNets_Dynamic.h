#ifndef BIONETS_DYNAMIC_H_INCLUDED
#define BIONETS_DYNAMIC_H_INCLUDED
#include "Analysis_BioNets_Linux.h"
//#include "Analysis_BioNets_Windows.h"
// Analysis of real biological networks.
int HammingDistance(int *v1,int *v2,int total){
    int sum=0;
    for(int ii=0;ii<total;ii++){
        sum+=(v1[ii]!=v2[ii]);}
    return sum;
}
int IsSame2Vector(int *v1,int *v2,int total){
    int logi=1;
    for(int ii=0;ii<total;ii++){
        if(v1[ii]!=v2[ii]){logi=0;break;}}
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
void FreeStatePrint(FILE *fp,int *ss,int total){
    fprintf(fp,"%1d",ss[0]);
    for(int ii=1;ii<total;ii++)fprintf(fp,"\t%1d",ss[ii]);
    fprintf(fp,"\n");
}

// sub-class of BioNets (logical)
class BioNetsDyn_Logic:public BioNets_Logical{
    public:
        int *StateVector;
        int *Candidate;
        int ***Loops;
        int *loop_num;
        void PrintGene(FILE *fp);
        void PrintState(FILE *fp);
        void FindAttractorLimitcycle(FILE *fp1,FILE *fp2,int TimeNum);
    private:
        int Update();
        void RandomVec();
        int Attractor();
        int LimitCyclePool_NewLoop(int index);
        int LimitCyclePool(int *RecordVec,int MaxLeng,int *GroupSize);
};
void BioNetsDyn_Logic::RandomVec(){
    for(int ii=0;ii<Gene_Network.total;ii++)StateVector[ii]=runif(mt)>0.5;
}
void BioNetsDyn_Logic::PrintGene(FILE *fp){
    fprintf(fp,"%s",GeneNames[0]);
    for(int ii=1;ii<Gene_Network.total;ii++)fprintf(fp,"\t%s",GeneNames[ii]);
    fprintf(fp,"\n");
}
void BioNetsDyn_Logic::PrintState(FILE *fp){
    fprintf(fp,"%1d",StateVector[0]);
    for(int ii=1;ii<Gene_Network.total;ii++)fprintf(fp,"\t%1d",StateVector[ii]);
    fprintf(fp,"\n");
}
int BioNetsDyn_Logic::Update(){
    int total=Gene_Network.total,*degs=Gene_Network.InDeg;
    int ii,jj=0,rand_id,tmps,tmp_int[total];
    for(ii=0;ii<total;ii++){
        if(Candidate[ii]){
            tmps=UpdateOneNode(StateVector,degs[ii],Parents[ii],Maps[ii]);
            if(tmps==StateVector[ii]){
                Candidate[ii]=0;}
            else {
                tmp_int[jj++]=ii;}}}
    if(jj>0){
        rand_id=(int)(runif(mt)*jj);
        rand_id=tmp_int[rand_id];// Random select one.
        StateVector[rand_id]=UpdateOneNode(StateVector,degs[rand_id],Parents[rand_id],Maps[rand_id]);
        struct Regulation *pt=Gene_Network.Network+rand_id;
        Candidate[rand_id]=0;
        while(pt->next!=NULL){
            Candidate[pt->next->code]=1;
            pt=pt->next;}}
    return jj;// "jj==0" means has reached an attractor.
}
int BioNetsDyn_Logic::Attractor(){
    int changes=1,step=0,MaxStep=1000;// MAX step:1000
    while(changes>0&&step<MaxStep){
        changes=BioNetsDyn_Logic::Update();
        step++;}
    return (step<MaxStep);// "Step<Max", succeed finding.
}
int BioNetsDyn_Logic::LimitCyclePool_NewLoop(int index){
    int total=Gene_Network.total,*degs=Gene_Network.InDeg;
    int number=1,**iptt,**NewLoop;
    int ii,jj=0,rand_id,tmps,logi_1,tmp_int[total],MaxPool=1024;
    struct Regulation *rpt2,*rpt,*source=(struct Regulation *)malloc(sizeof(struct Regulation));
    source->code=0;source->prev=source->next=NULL;
    NewLoop=(int**)malloc(sizeof(int*));
    NewLoop[0]=(int*)malloc(total*sizeof(int));
    memcpy(NewLoop[0],StateVector,total*sizeof(int));
    for(ii=0;ii<total;ii++){
        if(degs[ii]){
            tmps=UpdateOneNode(StateVector,degs[ii],Parents[ii],Maps[ii]);
            if(tmps!=StateVector[ii]){// Add candiate SS.
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
                number++;}}}
    while(source->next!=NULL&&number<MaxPool){
        memcpy(StateVector,NewLoop[source->next->code],total*sizeof(int));
        rpt=source->next;source->next=source->next->next;free(rpt);
        memcpy(tmp_int,StateVector,total*sizeof(int));
        for(ii=0;ii<total;ii++){
            if(degs[ii]&&number<MaxPool){
                tmps=UpdateOneNode(StateVector,degs[ii],Parents[ii],Maps[ii]);
                if(tmps!=StateVector[ii]){// Add candiate SS.
                    tmp_int[ii]=tmps;
                    logi_1=1;
                    for(jj=0;jj<number;jj++){
                        if(IsSame2Vector(NewLoop[jj],tmp_int,total)){logi_1=0;break;}}
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
                tmp_int[ii]=StateVector[ii];}}}
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
        number=-number;// Minus denotes a large loop.
    }
    return number;
}
int BioNetsDyn_Logic::LimitCyclePool(int *RecordVec,int MaxLeng,int *GroupSize){
    int ii,jj,dist,changes,logi,linshi;
    int ***ipttt,*ipt;int which_group;
    if(loop_num[0]==0){// The first loop-pool, small or large should be recorded.
        Loops=(int***)malloc(sizeof(int**));
        ipt=loop_num;
        loop_num=(int*)malloc(2*sizeof(int));
        memcpy(loop_num,ipt,sizeof(int));
        loop_num[1]=LimitCyclePool_NewLoop(loop_num[0]);
        loop_num[0]++;
        //printf("[%d:]\n",loop_num[0]);for(ii=0;ii<loop_num[0];ii++){printf("(%d)",loop_num[ii+1]);}printf("\n");
        /*printf("zhuyi!!!");getchar();
        for(int aa=0;aa<loop_num[1];aa++){for(int bb=0;bb<Gene_Network.total;bb++){fprintf(fpt,"%d",Loops[0][aa][bb]);}fprintf(fpt,"\n");}
        printf("[%p,%p]\n",Loops,Loops[0]);getchar();fclose(fpt);*/
        which_group=0;
        *GroupSize=loop_num[which_group+1];
    }
    else {
        logi=1;
        for(ii=0;ii<loop_num[0];ii++){
            if(logi){
                for(jj=0;jj<loop_num[ii+1];jj++){
                    if(IsSame2Vector(Loops[ii][jj],RecordVec,Gene_Network.total)){
                        which_group=ii;*GroupSize=loop_num[ii+1];
                        logi=0;break;}}}}
        /*printf("[%d:]",loop_num[0]);
        for(ii=0;ii<loop_num[0];ii++){
            printf("[%d]",loop_num[ii+1]);
        }printf("\n");*/
        if(logi==1){// New complex loop.
            ipttt=Loops;//printf("Has a new pool.\n");
            Loops=(int***)malloc((loop_num[0]+1)*sizeof(int**));
            memcpy(Loops,ipttt,sizeof(int**)*loop_num[0]);free(ipttt);
            //printf("[%d:]\n",loop_num[0]);
            linshi=BioNetsDyn_Logic::LimitCyclePool_NewLoop(loop_num[0]);
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
                    loop_num[0]++;
                }
                else {//Has existed!!
                    ipttt=Loops;//printf("Has a new pool.\n");
                    Loops=(int***)malloc((loop_num[0])*sizeof(int**));
                    memcpy(Loops,ipttt,sizeof(int**)*loop_num[0]);free(ipttt);
                    which_group=ii;}
                *GroupSize=-1024;}
            else {// Accept! A small pool.
                ipt=loop_num;
                loop_num=(int*)malloc((ipt[0]+2)*sizeof(int));
                memcpy(loop_num,ipt,(ipt[0]+1)*(sizeof(int)));free(ipt);
                loop_num[loop_num[0]+1]=linshi;
                which_group=loop_num[0];
                loop_num[0]++;
                *GroupSize=loop_num[loop_num[0]];}
            
            //printf("Has a new pool(okay,%p,%d).\n",Loops,loop_num[0]);
            //printf("[%d:]\n",loop_num[0]);
            //for(ii=0;ii<loop_num[0];ii++){printf("(%d)",loop_num[ii+1]);}getchar();
            }
        else ;}
    return which_group;// Regradless of repeated, record them.
}
void BioNetsDyn_Logic::FindAttractorLimitcycle(FILE *fp1,FILE *fp2,int TimeNum){
    int total=Gene_Network.total,Flag[total],Times;
    int loop_type,size[]={1};
    StateVector=(int*)malloc(total*sizeof(int));
    Candidate=(int*)malloc(total*sizeof(int));
    Loops=NULL;loop_num=(int *)malloc(sizeof(int));loop_num[0]=0;// Save each loop number.
    BioNetsDyn_Logic::PrintGene(fp1);
    fprintf(fp2,"period\n");
    if(total<24){// If system size is less than 24.
        Times=1<<total;
        for(int ii=0;ii<Times;ii++){
            memcpy(Candidate,Gene_Network.InDeg,total*sizeof(int));
            //Dec2Binary(StateVector,total,ii);
            BioNetsDyn_Logic::RandomVec();
            if(BioNetsDyn_Logic::Attractor()){// Is attractor, no need to saved.
                BioNetsDyn_Logic::PrintState(fp1);fprintf(fp2,"1\n");}
            else {//printf("Maybe a LC.");getchar();
                memcpy(Flag,StateVector,total*sizeof(int));
                loop_type=BioNetsDyn_Logic::LimitCyclePool(Flag,1000,size);
                fprintf(fp2,"%d\n",*size);// fprintf(fp2,"{%d,%d}\n",loop_type,*size);
                if(*size!=-1024)FreeStatePrint(fp1,Loops[loop_type][0],total);
                else BioNetsDyn_Logic::PrintState(fp1);}}}
    else {
        Times=TimeNum;
        for(int ii=0;ii<Times;ii++){
            memcpy(Candidate,Gene_Network.InDeg,total*sizeof(int));
            BioNetsDyn_Logic::RandomVec();
            if(BioNetsDyn_Logic::Attractor()){// Is attractor.
                BioNetsDyn_Logic::PrintState(fp1);fprintf(fp2,"1\n");}
            else {
                memcpy(Flag,StateVector,total*sizeof(int));
                loop_type=BioNetsDyn_Logic::LimitCyclePool(Flag,1000,size);
                fprintf(fp2,"%d\n",*size);// fprintf(fp2,"{%d,%d}\n",loop_type,*size);
                if(*size!=-1024)FreeStatePrint(fp1,Loops[loop_type][0],total);
                else BioNetsDyn_Logic::PrintState(fp1);}}}
    for(int ii=0;ii<loop_num[0];ii++){
        for(int jj=0;jj<loop_num[ii+1];jj++){
            free(Loops[ii][jj]);}
        free(Loops[ii]);}
    free(StateVector);free(Candidate);free(Loops);
    free(loop_num);
}
void exe_SearchAttractor(char *Dir,char *GeneNetCode,char *Output_Dir,char *ExpCodeNum,int Repeat){
    BioNetsDyn_Logic a_bionet;
    char target[64],ss[256],tt[512];
    FILE *fp_stab,*fp_step;
    std::string FilePathName;
    FilePathName.assign(Dir).append("//..//").append(Output_Dir).append("//").append(GeneNetCode);
    FilePathName.append("_").append(ExpCodeNum);
    strcpy(tt,FilePathName.c_str());
    fp_stab=fopen(strcat(tt,"_stab.txt"),"w");
    strcpy(tt,FilePathName.c_str());
    fp_step=fopen(strcat(tt,"_step.txt"),"w");
    FilePathName.erase(FilePathName.begin(),FilePathName.end());
    strcpy(tt,FilePathName.assign(Dir).append("//").append(GeneNetCode).c_str());
    a_bionet.GetBioNetwork(tt);
    a_bionet.FindAttractorLimitcycle(fp_stab,fp_step,Repeat);
    a_bionet.DeleteBioNetwork();
    fclose(fp_stab);fclose(fp_step);
}
#endif // BIONETS_DYNAMIC_H_INCLUDED