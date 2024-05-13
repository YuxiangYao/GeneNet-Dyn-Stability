/*
    This HEADER file is for Windows system. 
    Copyright (c) 2023-2024, YAO Yuxiang
    date: 2024-01-16
    version: 1.0
*/
#ifndef ANALYSIS_BIONETS_LINUX_H_INCLUDED
#define ANALYSIS_BIONETS_LINUX_H_INCLUDED
#include <sys/io.h> // Make sure your machine is X86 architecture, otherwise replace it by appropriate [file].h
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include "Boolean_Function.h"
#include "Network_GraphTheory.h"
#include <z3++.h>
// Remove comma from strings.
void Comma_0(char *ss,char sym){
    char *cpt=ss;
    while(1){
        if(*cpt==sym){
            *cpt='\0';break;}
        if(*cpt=='\0')break;
        cpt++;}
}
// Transfer logical to threshold-based formats.
int Z3_Logic_Threshold(int *maptab,int indeg,int length,int PlusMinus,FILE *fp){
    int ii,jj,res=0;
    z3::context TT_T;
    std::vector <z3::expr> vars;
    for(ii=0;ii<indeg;ii++){// Set variable number of VAR slots.
        vars.push_back(TT_T.int_const(("a_" + std::to_string(ii)).c_str()));}
    // Build solver.
    z3::expr theta=TT_T.int_const("theta");
    z3::expr exprss=theta;
    z3::solver solve_it(TT_T);
    for(ii=0;ii<length;ii++){
        exprss=theta;
        for(jj=0;jj<indeg;jj++){
            if((ii&(1<<jj))){// This bit is one (effective).
                exprss=exprss+vars[jj];}
            else if(PlusMinus)exprss=exprss-vars[jj];}
        if(maptab[ii]>0){
            solve_it.add(exprss>=0);}// The equation is .
        else {
            solve_it.add(exprss<=0);}
        for(jj=0;jj<indeg;jj++){// a_{xy}!=0
            solve_it.add(vars[jj]!=0);}}
    if(solve_it.check()==z3::sat){
        res=1;
        if(NULL!=fp){
            z3::model Models = solve_it.get_model();
            for(ii=0;ii<indeg;ii++){
                fprintf(fp,"%d,",Models.eval(vars[ii]).get_numeral_int());}// Each weighted links.
            fprintf(fp,"%d\n",Models.eval(theta).get_numeral_int());}}// theta: threshold.
    else {
        if(NULL!=fp){
            fprintf(fp,"Un-SAT [%d]\n",indeg);}}
    return res;
}
// Core class of biological networks.
class BioNets{
    public:
        char NetworkName[128];
        char **GeneNames; // Factor names.
        int **Maps; // All mapping table.
        int **Parents;// Parents slot.
        NetworkSys Gene_Network;// Structure of directed links.
        int ReturnGeneCode(char *item_name);// Return gene's code.
        void ShowBioNetwork();// Show the detail structure of bio-network.
        void DeleteBioNetwork();// Delete the bio-network and its slots.
        /* Analyzing functions */
        int TotalEdge();
        double TotalBias();
        double AverageSensitivity();// Return one bio-network's average sensitivity.
        void BoolFunFractionalAnalysis(FILE *fp);// Return fractions of Boolean Function's each type.
};
int BioNets::ReturnGeneCode(char *item_name){
    int code=0;
    while(1){
        if(strcmp(GeneNames[code],item_name)==0)break;
        code++;}
    return code;
}
void BioNets::ShowBioNetwork(){
    struct Regulation *npt,*net=Gene_Network.Network;
    int ii,total=Gene_Network.total;
    for(ii=0;ii<total;ii++,net++){
        npt=net;
        while(npt->prev!=NULL){
            printf("%s->",GeneNames[npt->prev->code]);npt=npt->prev;}
        npt=net;
        printf("[%s]",GeneNames[npt->code]);
        while(npt->next!=NULL){
            printf("->%s",GeneNames[npt->next->code]);npt=npt->next;}
        printf("\n");}
}
void BioNets::DeleteBioNetwork(){
    for(int ii=0;ii<Gene_Network.total;ii++){
        free(Maps[ii]);free(Parents[ii]);free(GeneNames[ii]);}
    free(GeneNames);GeneNames=NULL;
    free(Maps);Maps=NULL;
    free(Parents);Parents=NULL;
    Gene_Network.Delete_Network();
}
int BioNets::TotalEdge(){
    int ii,nn=Gene_Network.total;
    int *in_deg=Gene_Network.InDeg;
    int sum=0;
    for(ii=0;ii<nn;ii++){
        sum+=in_deg[ii];}
    return sum;
}
double BioNets::TotalBias(){
    int ii,jj,lens,nn=Gene_Network.total;
    int *in_deg=Gene_Network.InDeg;
    double sum=0,bias=0;
    for(ii=0;ii<nn;ii++){
        if(in_deg[ii]>0){
            lens=1<<in_deg[ii];
            sum+=(double)lens;
            for(jj=0;jj<lens;jj++){
                bias+=(double)(Maps[ii][jj]>0);}}}
    return bias/sum;
}
double BioNets::AverageSensitivity(){
    int ii,effective=0,nn=Gene_Network.total;
    int *in_deg=Gene_Network.InDeg;
    double sum=0,rrr;
    boolMap tmp_bool;
    for(ii=0;ii<nn;ii++){
        if(in_deg[ii]>1){// No-input nodes are ignored. 
            effective++;
            tmp_bool.k=in_deg[ii];
            tmp_bool.length=1<<in_deg[ii];
            tmp_bool.maptab=Maps[ii];
            rrr=tmp_bool.Sensitivity();//printf("[%s,%d]: %.4f\n",GeneNames[ii],in_deg[ii],rrr);
            sum+=rrr;}}
    return sum/effective;
}
void BioNets::BoolFunFractionalAnalysis(FILE *fp){
    int ii,RBF,tmp_if,nn=Gene_Network.total,*in_deg=Gene_Network.InDeg;
    boolMap_Cana tmp_c;    
    boolMap_Domin tmp_d;
    boolMap_Effec tmp_e;
    boolMap_Post tmp_p;
    boolMap_Unate tmp_u;
    boolMap_Thres tmp_t;
    for(ii=0;ii<nn;ii++){
        if(in_deg[ii]>1){// No-input nodes are ignored.
            RBF=0;
            fprintf(fp,"%s,%s,%d,",NetworkName,GeneNames[ii],in_deg[ii]);
            printf("%s,%s,%d\n",NetworkName,GeneNames[ii],in_deg[ii]);
            tmp_c.k=tmp_d.k=tmp_e.k=tmp_p.k=tmp_u.k=tmp_t.k=in_deg[ii];
            tmp_c.length=tmp_d.length=tmp_e.length=1<<in_deg[ii];
            tmp_p.length=tmp_u.length=tmp_t.length=1<<in_deg[ii];
            // IF cana?
            tmp_c.maptab=Maps[ii];tmp_if=tmp_c.is_Cana();RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            // IF dominate?
            tmp_d.Initial();tmp_d.In_Enumerate();tmp_d.Domin_Pre_Major();
            tmp_d.maptab=Maps[ii];tmp_if=tmp_d.is_Domin();RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            tmp_d.Reset();
            // IF effective?
            tmp_e.maptab=Maps[ii];tmp_if=tmp_e.is_Effec();fprintf(fp,"%d,",tmp_if);
            // IF post-2?
            tmp_p.maptab=Maps[ii];tmp_if=tmp_p.is_Post();RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            // IF monotanic? (Pure, Mixed)
            tmp_u.maptab=Maps[ii];tmp_if=tmp_u.is_Unate();RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            tmp_if=tmp_u.is_Mixed_Unate();RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            // IF threshold? (Standard, ) // St pn_1,theta=0,
            tmp_t.maptab=Maps[ii];tmp_if=tmp_t.is_Thres();RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            tmp_if=Z3_Logic_Threshold(Maps[ii],tmp_t.k,tmp_t.length,0,NULL);// 0/1
            RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            tmp_if=Z3_Logic_Threshold(Maps[ii],tmp_t.k,tmp_t.length,1,NULL);// +/-
            RBF+=tmp_if;fprintf(fp,"%d,",tmp_if);
            // Random.
            fprintf(fp,"%d\n",0==RBF);}}
}

// Sub class: logical type BioNet.
class BioNets_Logical:public BioNets{
    public:
        void GetBioNetwork(char *Dir); // Load a logical bio-network (in a folder).
        void Logic2Threshold(FILE *fp);// Obtain each bio-network's threshold-based formats.
    private:
        int ExternalComps(char *Dir_File,char **GeneSet); // Deal with external components.
        void LoadTruthTable(char *filename,char *Dir); // Load single genetic item with truth table.
};
int BioNets_Logical::ExternalComps(char *Dir_File,char **GeneSet){
    int sum=0;char target[64];target[0]='\0';
    FILE *fp=fopen(Dir_File,"r");
    while(!feof(fp)){
        fscanf(fp,"%s",target);
        if(target[0]=='\0')break;
        if(NULL!=GeneSet){
            GeneSet[sum]=(char*)malloc(64*sizeof(char));
            strcpy(GeneSet[sum],target);}
        target[0]='\0';
        sum++;}
    fclose(fp);
    return sum;
}
void BioNets_Logical::GetBioNetwork(char *Dir){
    int ii,SysSize=0,externs=0;
    char external_comp[128];//single_name[128];
    std::string FilePath;
    std::vector<std::string> FileList;
    struct dirent *FileInfo;// Structor of file information.
    DIR *dir_ptr;
    dir_ptr=opendir(FilePath.assign(Dir).append("//").c_str());
    if(dir_ptr!=NULL){
        while((FileInfo=readdir(dir_ptr))!=NULL){
            if(strcmp(FileInfo->d_name, ".")!=0&&strcmp(FileInfo->d_name, "..")!=0){
                if(strcmp(FileInfo->d_name,"external_components.ALL.txt")!=0){
                    FileList.push_back(FilePath.assign(FileInfo->d_name));
                    SysSize++;}
                else {// External components filename.
                    strcpy(external_comp,
                        FilePath.assign(Dir).append("//").append(FileInfo->d_name).c_str());}}}
        sort(FileList.begin(),FileList.end());// Must be sorted!!
        closedir(dir_ptr);
        // Set slot of network.
        GeneNames=NULL;// Set as NULL for only counting external components.
        externs=BioNets_Logical::ExternalComps(external_comp,GeneNames);
        Gene_Network.total=SysSize+externs;// All components: core + extern
        // Record genenames.
        GeneNames=(char**)malloc(Gene_Network.total*sizeof(char*));
        Maps=(int**)malloc(Gene_Network.total*sizeof(int*));
        Parents=(int**)malloc(Gene_Network.total*sizeof(int*));
        Gene_Network.Initial_Network();// Initial the network.
        externs=BioNets_Logical::ExternalComps(external_comp,GeneNames+SysSize);// Exntern in last.
        for(ii=0;ii<externs;ii++){// External components have not inputs.
            *(Parents+SysSize+ii)=NULL; *(Maps+SysSize+ii)=NULL;}
        // Remodified genes' name.
        for(ii=0;ii<SysSize;ii++){
            GeneNames[ii]=(char*)malloc(64*sizeof(char));
            strcpy(GeneNames[ii],FileList[ii].c_str());
            Comma_0(GeneNames[ii],'.');}
        for(ii=0;ii<SysSize;ii++){
            BioNets_Logical::LoadTruthTable(GeneNames[ii],Dir);}}
    else {
        printf("Error in reading files.\n");getchar();}
}
void BioNets_Logical::LoadTruthTable(char *filename,char *Dir){
    FILE *fp;
    char tmp_name[256],target[256];
    strcpy(tmp_name,Dir);strcat(tmp_name,"//");
    strcat(tmp_name,filename);strcat(tmp_name,".csv");
    int ii,pars=0,length,*tmp_map;
    int parent,Child;
    if((fp=fopen(tmp_name,"r"))==NULL){
        printf("The %s not exist!\n",tmp_name);getchar();}
    else {
        Child=BioNets::ReturnGeneCode(filename);// Child's code number.
        pars=1;
        while(!feof(fp)&&pars){
            fgets(target,sizeof(target),fp);// Read first line.
            ii=0;length=0;
            while(target[ii]!='\0'){
                if(target[ii]==',')length++;
                ii++;}
            pars=0;}
        rewind(fp);
        Parents[Child]=(int*)malloc(length*sizeof(int));
        while(!feof(fp)){
            while(1){
                fscanf(fp,"%s",target);
                Comma_0(target,',');
                if(target[0]=='1'||target[0]=='0')break;
                if(pars==length)break;
                parent=BioNets::ReturnGeneCode(target);
                Parents[Child][length-pars-1]=parent;
                BuildRegulationship(Gene_Network.Network,Child,parent);
                (Gene_Network.OtDeg)[parent]++;
                pars++;}}
        (Gene_Network.InDeg)[Child]=pars;
        ii=0;
        rewind(fp);
        fgets(target,sizeof(target),fp);// Ignore the title line.
        tmp_map=Maps[Child]=(int*)malloc(sizeof(int)<<pars);
        pars++;// Include target node.
        while(!feof(fp)){
            ii++;
            fscanf(fp,"%s",target);
            if((ii%pars)==0){
                *tmp_map++=atoi(target);}}}
    fclose(fp);
}
void BioNets_Logical::Logic2Threshold(FILE *fp){
    int ii,nn=Gene_Network.total,tmp;
    int *in_deg=Gene_Network.InDeg;
    int *ot_deg=Gene_Network.OtDeg;
    for(ii=0;ii<nn;ii++){
        if(in_deg[ii]>0&&ot_deg[ii]>0){// No-input nodes are ignored. 
            fprintf(fp,"%s: ",GeneNames[ii]);printf("%s >>",GeneNames[ii]);
            tmp=Z3_Logic_Threshold(Maps[ii],in_deg[ii],(1<<in_deg[ii]),0,fp);}}
    printf("\n");
}

// Sub_Class: threshold-based models.
class BioNets_Threshold:public BioNets{
    public:
        void LoadAllNodes(char *dir_filename);// Threshold based BF network analysis.
        void LoadThresholdEdges(char *dir_filename);// Load the threshold based network.
    private:
        void SetNewNode(int **LinkType,int *InDeg,int target,int type);// Set new topology.
        void AttachNode(int **LinkType,int *InDeg,int target,int type);// Add at existed topology.
};
void BioNets_Threshold::LoadAllNodes(char *dir_filename){
    FILE *fp;
    int logi=0,ii,*ipt;
    char tmp_name[100],node[64],num[20],*cpt;
    strcpy(tmp_name,dir_filename);
    if((fp=fopen(strcat(tmp_name,".ids"),"r"))==NULL){
        printf("The %s not exist!\n",tmp_name);getchar();}
    while(!feof(fp)){// Count number of nodes.
        cpt=fgets(node,sizeof(node),fp);
        if(cpt!=NULL)logi++;}
    // Initialize some slots.
    Gene_Network.total=logi;// All components.
    GeneNames=(char**)malloc(logi*sizeof(char*));
    Maps=(int**)malloc(logi*sizeof(int*));
    Parents=(int**)malloc(logi*sizeof(int*));
    Gene_Network.Initial_Network();
    rewind(fp);
    for(ii=0;ii<logi;ii++){
        Parents[ii]=NULL;// Unuse for threshold-based.
        GeneNames[ii]=(char*)malloc(64*sizeof(char));// Max char-input: 64
        fscanf(fp,"%s",node);fscanf(fp,"%s",num);
        strcpy(GeneNames[ii],node);}
    fclose(fp);
}
void BioNets_Threshold::SetNewNode(int **LinkType,int *InDeg,int target,int type){
    int *ipt;
    InDeg[target]++;
    ipt=(int*)malloc(sizeof(int));
    *ipt=type;
    LinkType[target]=ipt;
}
void BioNets_Threshold::AttachNode(int **LinkType,int *InDeg,int target,int type){
    int *ipt,INT_BIT=sizeof(int),degree=InDeg[target],length=INT_BIT*degree;
    ipt=(int*)malloc(length+INT_BIT);
    memcpy(ipt,LinkType[target],length);
    ipt[degree]=type;
    free(LinkType[target]);
    LinkType[target]=ipt;
    InDeg[target]++;
}
void BioNets_Threshold::LoadThresholdEdges(char *dir_filename){
    FILE *fp;
    char tmp_name[128],target[128],source[128],type[16];
    strcpy(tmp_name,dir_filename);
    strcat(tmp_name,".topo");
    int ii,jj,kk,sum,length,parent,child,*tmp_map,total=Gene_Network.total;
    int *ind=Gene_Network.InDeg,*otd=Gene_Network.OtDeg;
    int **topo_link;
    if((fp=fopen(tmp_name,"r"))==NULL){
        printf("The %s not exist!\n",tmp_name);getchar();}
    else {
        topo_link=(int**)malloc(total*sizeof(int*));
        for(ii=0;ii<total;ii++){
            topo_link[ii]=NULL;}
        fgets(source,sizeof(source),fp);// Ignore title
        while(!feof(fp)){
            kk=fscanf(fp,"%s",source);
            fscanf(fp,"%s",target);
            fscanf(fp,"%s",type);
            if(atoi(type)<2){sum=1;}else {sum=-1;}
            if(kk>0){// Valid read.
                child=BioNets::ReturnGeneCode(target);
                parent=BioNets::ReturnGeneCode(source);
                BuildRegulationship(Gene_Network.Network,child,parent);// For strucutre gene.network.
                otd[parent]++;
                if(ind[child]>0){// The indeg has been updated.
                    AttachNode(topo_link,ind,child,sum);}
                else {
                    SetNewNode(topo_link,ind,child,sum);}}}
        // Now transform +/- topology into truth map.
        for(ii=0;ii<total;ii++){
            Maps[ii]=NULL;
            if(ind[ii]>0){
                length=(1<<ind[ii]);// WARNING! In-Deg=0 is ignored, NULL;
                tmp_map=(int*)malloc(length*sizeof(int));
                for(jj=0;jj<length;jj++){
                    sum=0;
                    for(kk=0;kk<ind[ii];kk++){// Turn 0/1 to +/-.
                        //sum+=topo_link[ii][kk]*((((jj&1<<kk)>0)<<1)-1);// +1/-1
                        sum+=topo_link[ii][kk]*((jj&1<<kk)>0);}// 0/1
                    tmp_map[jj]=sum>0;}
                Maps[ii]=tmp_map;free(topo_link[ii]);}}
        free(topo_link);}
    fclose(fp);
}

// General properties of all bio-networks. (Av. Sensitivity)
void exe_AvSensitivity_Logic(char *Dir,char *FileList,char *Output_File){
    BioNets_Logical a_bionet;
    char target[64],ss[128];
    strcpy(ss,Dir);strcat(ss,FileList);
    FILE *fp_list=fopen(ss,"r");
    strcpy(ss,Dir);strcat(ss,Output_File);
    FILE *fp_results=fopen(ss,"w");
    while(!feof(fp_list)){
        fscanf(fp_list,"%s",target);
        if(feof(fp_list))break;
        strcpy(ss,Dir);strcat(ss,target);
        a_bionet.GetBioNetwork(ss);//a_bionet.ShowBioNetwork();getchar();
        fprintf(fp_results,"Net_%05d\t%d\t",atoi(target),a_bionet.Gene_Network.total);
        fprintf(fp_results,"%d\t%.4f\t",a_bionet.TotalEdge(),a_bionet.TotalBias());
        fprintf(fp_results,"%.4f\n",a_bionet.AverageSensitivity());
        a_bionet.DeleteBioNetwork();}
    fclose(fp_list);fclose(fp_results);
}
void exe_AvSensitivity_Thres(char *Dir,char *FileList,char *Output_File){
    BioNets_Threshold a_bionet;
    char target[64],ss[128];
    strcpy(ss,Dir);strcat(ss,FileList);
    FILE *fp_list=fopen(ss,"r");
    strcpy(ss,Dir);strcat(ss,Output_File);
    FILE *fp_results=fopen(ss,"w");
    while(!feof(fp_list)){
        fscanf(fp_list,"%s",target);
        strcpy(ss,Dir);
        strcat(ss,target);
        strcpy(a_bionet.NetworkName,"N_");
        strcat(a_bionet.NetworkName,target);
        a_bionet.LoadAllNodes(ss);
        a_bionet.LoadThresholdEdges(ss);
        fprintf(fp_results,"Net_%s\t%d\t",a_bionet.NetworkName,a_bionet.Gene_Network.total);
        fprintf(fp_results,"%d\t%.4f\t",a_bionet.TotalEdge(),a_bionet.TotalBias());
        fprintf(fp_results,"%.4f\n",a_bionet.AverageSensitivity());
        a_bionet.DeleteBioNetwork();}
    fclose(fp_list);
    fclose(fp_results);
}
// Obtain each bio-network's function fractional components.
void exe_FunFracComp_Logic(char *Dir,char *FileList,char *Output_File){
    BioNets_Logical a_bionet;
    char target[64],ss[128];
    strcpy(ss,Dir);strcat(ss,FileList);
    FILE *fp_list=fopen(ss,"r");
    strcpy(ss,Dir);strcat(ss,Output_File);
    FILE *fp_results=fopen(ss,"w");
    while(!feof(fp_list)){
        fscanf(fp_list,"%s",target);
        if(feof(fp_list))break;
        strcpy(ss,Dir);
        strcat(ss,target);
        strcpy(a_bionet.NetworkName,"Net_");
        strcat(a_bionet.NetworkName,target);
        printf("%s\n",a_bionet.NetworkName);
        a_bionet.GetBioNetwork(ss);
        a_bionet.BoolFunFractionalAnalysis(fp_results);
        a_bionet.DeleteBioNetwork();}
    fclose(fp_list);
    fclose(fp_results);
}
void exe_FunFracComp_Thres(char *Dir,char *FileList,char *Output_File){
    BioNets_Threshold a_bionet;
    char target[64],ss[128];
    strcpy(ss,Dir);strcat(ss,FileList);
    FILE *fp_list=fopen(ss,"r");
    strcpy(ss,Dir);strcat(ss,Output_File);
    FILE *fp_results=fopen(ss,"w");
    int ii,jj;
    while(!feof(fp_list)){
        fscanf(fp_list,"%s",target);
        strcpy(ss,Dir);
        strcat(ss,target);
        strcpy(a_bionet.NetworkName,"N_");
        strcat(a_bionet.NetworkName,target);
        a_bionet.LoadAllNodes(ss);
        a_bionet.LoadThresholdEdges(ss);
        a_bionet.BoolFunFractionalAnalysis(fp_results);
        a_bionet.DeleteBioNetwork();}
    fclose(fp_list);fclose(fp_results);
}
// Transfer the logical to the threshold.
void exe_Logic_to_Threshold(char *Dir,char *FileList,char *Output_File){
    BioNets_Logical a_bionet;
    char target[64],ss[128];
    strcpy(ss,Dir);strcat(ss,FileList);
    FILE *fp_list=fopen(ss,"r");
    strcpy(ss,Dir);strcat(ss,Output_File);
    FILE *fp_results=fopen(ss,"w");
    while(!feof(fp_list)){
        fscanf(fp_list,"%s",target);
        if(target[0]=='\0')break;
        strcpy(ss,Dir);
        strcat(ss,target);
        a_bionet.GetBioNetwork(ss);
        fprintf(fp_results,"Net_%s\n",target);
        printf("#%s\n",target);
        a_bionet.Logic2Threshold(fp_results);
        fprintf(fp_results,"\n");
        a_bionet.DeleteBioNetwork();}
    fclose(fp_list);fclose(fp_results);
}
// Random testing abundence of threshold.
void exe_RandomTesting(){
    int ii,k,length,*maps;FILE *fp=NULL;
    int ccc=0,nums[16],user[16];
    memset(nums,0,16*sizeof(int));
    memset(user,0,16*sizeof(int));
    while(ccc<500){
        k=(int)(runif(mt)*7)+2;
        length=1<<k;
        maps=(int*)malloc(length*sizeof(int));
        for(ii=0;ii<length;ii++){
            maps[ii]=r_binary(mt);}//std::cout << "[" <<k<<"] ";
        ii=Z3_Logic_Threshold(maps,k,length,0,fp);
        free(maps);maps=NULL;
        user[k]++;
        if(ii>0){//std::cout<< "Yes" <<std::endl;
            nums[k]++;}
        else ;//std::cout<< "No" <<std::endl;
        ccc++;}
    for(ii=0;ii<16;ii++){
        printf("[%d]: %d/%d\n",ii,nums[ii],user[ii]);
    }getchar();
}
#endif // ANALYSIS_BIONETS_LINUX_H_INCLUDED
