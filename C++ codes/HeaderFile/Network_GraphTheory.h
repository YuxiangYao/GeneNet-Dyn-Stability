/*
  Copyright (c) 2022-2024, YAO Yuxiang
    date: 2023-04-18
    version: 2.0
  Network generator:
*/
#ifndef NETWORK_GRAPHTHEORY_H_INCLUDED
#define NETWORK_GRAPHTHEORY_H_INCLUDED
#include "Common_Used.h"
struct Regulation{// Structure of regulations.
    int code;// Node's code.
    struct Regulation *prev;// The previous node.
    struct Regulation *next;// The next node.
};
void InsertParentsNode(struct Regulation *p, int parent){
    struct Regulation *tmp=(struct Regulation *)malloc(sizeof(struct Regulation));
    tmp->code=parent;
    if(p->prev==NULL){
        p->prev=tmp;tmp->prev=NULL;tmp->next=p;}
    else {
        p->prev->next=tmp;tmp->prev=p->prev;tmp->next=p;p->prev=tmp;}
}
void InsertChildsNode(struct Regulation *p,int child){
    struct Regulation *tmp=(struct Regulation *)malloc(sizeof(struct Regulation));
    tmp->code=child;
    if(p->next==NULL){
        p->next=tmp;tmp->prev=p;tmp->next=NULL;}
    else {
        p->next->prev=tmp;tmp->next=p->next;p->next=tmp;tmp->prev=p;}
}
void BuildRegulationship(struct Regulation *net,int child,int parent){
    InsertParentsNode(&net[child],parent);
    InsertChildsNode(&net[parent],child);
}
// Build non-directional edge (Only use child, not contain the parent-type)
void BuildNeighborEdge(struct Regulation *net,int node1,int node2){
    InsertChildsNode(&net[node1],node2);
    InsertChildsNode(&net[node2],node1);
}
// Remove one specified child from parent node.
int DeleteChilds(struct Regulation *p,int child){
    struct Regulation *npt=p,*Dels;
    int GotIt=0;
    while(npt->next!=NULL){
        if(npt->next->code==child){// Find where target node
            GotIt=1;
            if(npt->next->next==NULL){
                free(npt->next);npt->next=NULL;}
            else {
                npt->next->next->prev=npt;Dels=npt->next;
                npt->next=npt->next->next;free(Dels);}
            break;}
        else npt=npt->next;}
    return GotIt;
}
// Remove one specified parent from child node.
int DeleteParents(struct Regulation *p,int parent){
    struct Regulation *npt=p,*Dels;int GotIt=0;
    while(npt->prev!=NULL){// Delete parent target in this son node
        if(npt->prev->code==parent){// Find where target vertex
            GotIt=1;
            if(npt->prev->prev==NULL){
                free(npt->prev);npt->prev=NULL;}
            else {
                npt->prev->prev->next=npt;Dels=npt->prev;
                npt->prev=npt->prev->prev;free(Dels);}
            break;}
        else npt=npt->prev;}
    return GotIt;
}
// Node has No.xx child/parent.
int ExistChild(struct Regulation *p,int child){
    struct Regulation *npt=p; int inhere=0;
    while(npt->next!=NULL){
        if(child==npt->next->code){inhere=1;break;}
        else npt=npt->next;}
    return inhere;
}
int ExistParent(struct Regulation *p,int parent){
    struct Regulation *npt=p; int inhere=0;
    while(npt->prev!=NULL){
        if(parent==npt->prev->code){inhere=1;break;}
        else npt=npt->prev;}
    return inhere;
}
// Build a mutual loop of all nodes.
void LinkMutualLoop(struct Regulation *Loops,int total){
    int ii;struct Regulation *rpt=Loops;
    for(ii=0;ii<total;ii++,rpt++)rpt->code=ii;
    Loops[0].next=Loops+1;
    Loops[0].prev=Loops+total-1;
    Loops[total-1].next=Loops;
    Loops[total-1].prev=Loops+total-2;
    for(ii=1;ii<total-1;ii++){
        Loops[ii].next=Loops+ii+1;
        Loops[ii].prev=Loops+ii-1;}
}
void LinkRomoveCurrent(struct Regulation *rpt){
    rpt->prev->next=rpt->next;
    rpt->next->prev=rpt->prev;
    rpt->next=rpt->prev=rpt;// Current node points itself!
}
// First-in first-out to queue of distance.
void ShortestFind_FIFO(struct Regulation *ParentNode,int Code,int logi){// For shortest path algorithm.
    struct Regulation *tmp=ParentNode;
    if(logi>0){// Insert candidate nodes at end.
        struct Regulation *Newer=(struct Regulation *)malloc(sizeof(struct Regulation));
        Newer->next=Newer->prev=NULL;
        Newer->code=Code;
        while(tmp->next!=NULL){// Follow the last member.
            tmp=tmp->next;}
        tmp->next=Newer;
        Newer->prev=tmp;}
    else {// Delete the queue's first.
        if(ParentNode->next->next==NULL){// Only one member.
            tmp=ParentNode->next;
            ParentNode->next=NULL;}
        else {
            tmp=tmp->next;
            ParentNode->next->next->prev=ParentNode;
            ParentNode->next=ParentNode->next->next;}
        free(tmp);}
}
// Auxiliary function for lattice undirected.
void Auxiliary_LT_U(int Length,int Height,int Type,int ids,int *Neighbor){// Note: ids only 3, 4, 6; total is a square number (even*even).
    int ii=ids/Length,jj=ids%Length,tmp1,tmp2,tmp3;
    if(Type==4){// Square
        Neighbor[0]=jj-1;if(Neighbor[0]<0)Neighbor[0]=Length-1;Neighbor[0]=Neighbor[0]+ii*Length;
        Neighbor[1]=jj+1;if(Neighbor[1]==Length)Neighbor[1]=0;Neighbor[1]=Neighbor[1]+ii*Length;
        Neighbor[2]=ii-1;if(Neighbor[2]<0)Neighbor[2]=Height-1;Neighbor[2]=Neighbor[2]*Length+jj;
        Neighbor[3]=ii+1;if(Neighbor[3]==Height)Neighbor[3]=0;Neighbor[3]=Neighbor[3]*Length+jj;}
    else if(Type==3){// Triangle
        Neighbor[0]=jj-1;if(Neighbor[0]<0)Neighbor[0]=Length-1;Neighbor[0]=Neighbor[0]+ii*Length;
        Neighbor[1]=jj+1;if(Neighbor[1]==Length)Neighbor[1]=0;Neighbor[1]=Neighbor[1]+ii*Length;
        if((ii+jj)%2==0){// even
            Neighbor[2]=ii+1;if(Neighbor[2]==Height)Neighbor[2]=0;Neighbor[2]=Neighbor[2]*Length+jj;}
        else {// odd
            Neighbor[2]=ii-1;if(Neighbor[2]<0)Neighbor[2]=Height-1;Neighbor[2]=Neighbor[2]*Length+jj;}}
    else if(Type==6){// Hexagon
        Neighbor[0]=jj-1;if(Neighbor[0]<0)Neighbor[0]=Length-1;Neighbor[0]=Neighbor[0]+ii*Length;
        Neighbor[1]=jj+1;if(Neighbor[1]==Length)Neighbor[1]=0;Neighbor[1]=Neighbor[1]+ii*Length;
        tmp1=ii-1;tmp2=ii+1;tmp3=jj-1;
        if(tmp1<0)tmp1=Height-1;if(tmp2==Height)tmp2=0;if(tmp3<0)tmp3=Length-1;
        Neighbor[2]=tmp1*Length+jj;
        Neighbor[3]=tmp1*Length+tmp3;
        Neighbor[4]=tmp2*Length+jj;
        Neighbor[5]=tmp2*Length+tmp3;}
}

// Define network class.
class NetworkSys {
    public:
        int total;// Size of network
        int *InDeg;
        int *OtDeg;
        struct Regulation *Network;// Slot for edge's relations.
        // Basic function:
        void Initial_Network();
        void Reset_Network();
        void Delete_Network();
        void Show_Network();
        void Count_InOt_Deg();
        void Show_Degree();
        void ShortestPath(int source);
        // Generate various types of network
        //void Net_ER_U1(double AvDeg);
        //void Net_ER_U2(double Prob);
        //void Net_ER_D(double AvDeg);
        //void Net_SF_U(int MaxDeg,int MinDeg,double Gamma);
        //void Net_BA_U(int Core);
        //void Net_RR_U(int K);
        //void Net_RR_D3(int K);
        //void Net_RR_D1(int K);
        //void Net_LT_U(int K);
        //void Net_NK_D(int K);
        //void Net_Fr_U(int K);
        //void Net_Cy_U(int K);
    private:
        int *Coding;// Code: 0,1,2,3,4,...
        int **Address;// Pointing address of [code] (Double of total+1).
        void Address_Reset();// Reset temporary address slot. 
};
// Initial one network.[V]
void NetworkSys::Initial_Network(){
    Network=(struct Regulation *)malloc(total*sizeof(struct Regulation));
    InDeg=(int *)malloc(4*total*sizeof(int));// InDeg+OtDeg+Coding+ReserveSlot
    OtDeg=InDeg+total;
    Coding=OtDeg+total;
    Address=(int **)malloc(2*(total+1)*sizeof(int *));
    int **iptt=Address+total+1;
    struct Regulation *npt=Network;
    for(int ii=0;ii<total;ii++,npt++){
        InDeg[ii]=OtDeg[ii]=0;
        Coding[ii]=ii;
        iptt[ii]=Coding+ii;
        npt->code=ii;npt->prev=npt->next=NULL;}
    iptt[total]=NULL;
}
// Reset one network (Only delete and free the relations).[V]
void NetworkSys::Reset_Network(){
    struct Regulation *del,*npt=Network;
    for(int ii=0;ii<total;ii++,npt++){
        InDeg[ii]=OtDeg[ii]=0;
        while(npt->prev!=NULL){// Delete parents
            del=npt->prev;npt->prev=npt->prev->prev;free(del);}
        while(npt->next!=NULL){// Delete children
            del=npt->next;npt->next=npt->next->next;free(del);}}
}
// Delete one network (Free all memory).[V]
void NetworkSys::Delete_Network(){
    NetworkSys::Reset_Network();// First reset.
    free(InDeg);InDeg=OtDeg=Coding=NULL;
    free(Network);Network=NULL;
    free(Address);Address=NULL;
}
// Show one network.[V]
void NetworkSys::Show_Network(){
    struct Regulation *npt,*net=Network;
    for(int ii=0;ii<total;ii++,net++){
        npt=net;
        while(npt->prev!=NULL){
            printf("%d->",npt->prev->code);npt=npt->prev;}
        npt=net;
        printf("[%d]",npt->code);
        while(npt->next!=NULL){
            printf("->%d",npt->next->code);npt=npt->next;}
        printf("\n");} 
}
// Show each in/out degree.[V]
void NetworkSys::Show_Degree(){
    int *ind=InDeg,*otd=OtDeg;
    for(int ii=0;ii<total;ii++,ind++,otd++){
        printf("[%d]  %d,%d\n",ii,*ind,*otd);}
}
// Count in/out degree of each vertex.[V]
void NetworkSys::Count_InOt_Deg(){
    struct Regulation *net=Network,*npt;
    for(int ii=0;ii<total;ii++,net++){
        npt=net;
        while(npt->next!=NULL){// Count "npt->prev" is same for directional net.
            InDeg[npt->next->code]++;
            OtDeg[ii]++;
            npt=npt->next;}}
}
// Distance table of one node to all (Dijkstra, single source).
void NetworkSys::ShortestPath(int source){// Priority queue, not save routes
    int tmpI,tmp_source,tmp_node;
    int *Distance=OtDeg+(total<<1);
    struct Regulation *pt=NULL,*head=(struct Regulation *)malloc(sizeof(struct Regulation));
    head->next=head->prev=NULL;head->code=-1;
    for(tmpI=0;tmpI<total;tmpI++){Distance[tmpI]=total<<1;}// Initialization
    Distance[source]=0;
    ShortestFind_FIFO(head,source,666);// Insert source node
    if(Network[source].next!=NULL){
        while(1){
            if(head->next==NULL)break;
            tmp_source=head->next->code;
            ShortestFind_FIFO(head,-888,-666);// Out of queue (2nd Argu useless)
            pt=&Network[tmp_source];
            while(pt->next!=NULL){
                tmp_node=pt->next->code;
                if(Distance[tmp_source]+1<Distance[tmp_node]){// A shorter path
                    Distance[tmp_node]=Distance[tmp_source]+1;
                    ShortestFind_FIFO(head,tmp_node,666);}// Insert candidate node
                pt=pt->next;}}}
    free(head);
}
// Rest address.[V]
void NetworkSys::Address_Reset(){
    memcpy(Address,Address+total+1,(total+1)*sizeof(int *));
}
#endif // NETWORK_GRAPHTHEORY_H_INCLUDED
