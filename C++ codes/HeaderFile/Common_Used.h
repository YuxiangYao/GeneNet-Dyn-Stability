#ifndef COMMON_USED_H_INCLUDED
#define COMMON_USED_H_INCLUDED
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h> 
#include <algorithm>
#include "Random_Setting.h"
// Unrepeatably choose K elements from N elements. 
void Choose_NK(int N,int K,int *Selected,int *Label,int flag){
    /*  Selected: record the select IDs 
        Label: labels for recording which node has been chosen 
        flag: Can select pointed node?   */
    int ii;
    if(K<N/2){// Small selection.
        int tmp;
        if(flag>=0)Label[flag]=1;
        for(ii=0;ii<K;ii++){
            do{
                tmp=(int)(N*runif(mt));
            }while(Label[tmp]);
            Selected[ii]=tmp;
            Label[tmp]=1;}
        if(flag>=0){// Recover the slot Label
            Label[flag]=0;}
        for(ii=0;ii<K;ii++){
            Label[Selected[ii]]=0;}}
    else {// Large selection.
        int *ipt,elements[N];
        ipt=elements;
        for(ii=0;ii<N;ii++,ipt++)*ipt=ii;
        std::shuffle(elements,elements+N,mt);
        memcpy(Selected,elements,sizeof(int)*K);}// Only choose former K.
}
// Calculate the average and standard deviation of vector "vars".
void Mean_SD(double *vars,int total,double *results){
    /*  vars: Variables to be solved.
        total: number of Variables.
        results: vectors with 2 elements [0]:Av., [1]:sd. */
    int tmpI;double sum=0;
    for(tmpI=0;tmpI<total;tmpI++)sum+=vars[tmpI];
    results[0]=sum/((double)total);sum=0;
    for(tmpI=0;tmpI<total;tmpI++){
        sum+=(vars[tmpI]-results[0])*(vars[tmpI]-results[0]);}
    results[1]=sqrt(sum/((double)total));
}
// Link all strings by given hyphen.
void StringPaste(char **ss,char *hyphen,int num){
    for(int ii=1;ii<num;ii++){// [NOTE] s0 must be sufficiently large.
        strcat(ss[0],hyphen); 
        strcat(ss[0],ss[ii]);}
}
#endif // COMMON_USED_H_INCLUDED
