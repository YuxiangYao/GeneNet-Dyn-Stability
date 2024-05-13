#ifndef BOOLEAN_FUNCTION_SENSITIVITY_H_INCLUDED
#define BOOLEAN_FUNCTION_SENSITIVITY_H_INCLUDED
#include "Boolean_Function.h"
// Class of sensitivity analysis (Boolean).
class B_Sensitivity_Analysis{
    public:
        char *ProjectName;
        int Loop;// Repeat times for analyzing
        FILE *fpss;// Output file's pointer.
        void Initialization(char *Type,char *Times,char *RandSeed);
        void Delete();// Delete, reset, and free momery.
        void Subtype(char *BoolFunType,char *Times,char *RandSeed,char *SubSubType,char *OtherArgu);
};
void B_Sensitivity_Analysis::Initialization(char *Type,char *Times,char *RandSeed){
    mt.seed(atoi(RandSeed));
    Loop=atoi(Times);
    ProjectName=(char *)malloc(128*sizeof(char));
    char *cc[8],tmp[2]={"_"},FileName[128]={"b_Sensi"};
    memcpy(ProjectName,FileName,128*sizeof(char));
    cc[0]=ProjectName;
    cc[1]=Type;
    cc[2]=Times;
    cc[3]=RandSeed;
    StringPaste(cc,tmp,4);
}
void B_Sensitivity_Analysis::Delete(){
    free(ProjectName);
    ProjectName=NULL;
}
void B_Sensitivity_Analysis::Subtype(char *BoolFunType,char *Times,char *RandSeed,char *SubSubType,char *OtherArgu){
    B_Sensitivity_Analysis::Initialization(BoolFunType,Times,RandSeed);
    char *cc[8],tmp[2];tmp[0]='_';cc[0]=ProjectName;cc[1]=SubSubType;cc[2]=OtherArgu;
    switch(BoolFunType[0]){// Configurate the parameters according to various types.
        case 'R':case 'E':break;case 'C':StringPaste(cc,tmp,3);break;
        case 'D':case 'P':case 'U':case 'T':StringPaste(cc,tmp,2);break;}
    fpss=fopen(strcat(ProjectName,".txt"),"w");
    int ii,jj,*Maps=NULL,Length=4;// First length is 2^2=4.
    double m_sen,sensty_s[Loop<<2],energy_s[Loop<<2],io[2];
    double *sensty_r,*energy_r,*sensty__,*energy__,*MajorBia;
    sensty_r=sensty_s+Loop;sensty__=sensty_r+Loop;MajorBia=sensty__+Loop;
    energy_r=energy_s+Loop;energy__=energy_r+Loop;
    boolMap BFs;        boolMap RBFs; 
    boolMap_Cana  CBFs; boolMap_Domin DBFs;
    boolMap_Effec EBFs; boolMap_Post  PBFs;
    boolMap_Thres TBFs; boolMap_Unate UBFs;
    for(ii=2;ii<4;ii++){// K gets from 2 to 10.
        Maps=(int*)malloc(Length*sizeof(int));
        m_sen=0.500*ii;
        BFs.k=ii;BFs.length=Length;BFs.maptab=Maps;BFs.Initial();
        switch(BoolFunType[0]){
            case 'C':CBFs.k=ii;  CBFs.length=Length;  CBFs.maptab=Maps;
                    CBFs.Deep=atoi(SubSubType);  CBFs.Fixed=atoi(OtherArgu);
                    CBFs.Cana_Initial();  CBFs.In_Enumerate();  break;
            case 'D':DBFs.k=ii;  DBFs.length=Length;  DBFs.Initial();  DBFs.maptab=Maps;
                    DBFs.Domin_Type=atoi(SubSubType);  DBFs.In_Enumerate();  DBFs.Domin_Pre_Major();  break;
            case 'E':EBFs.k=ii;  EBFs.length=Length;  EBFs.Initial();  EBFs.maptab=Maps;  break;
            case 'P':PBFs.k=ii;  PBFs.length=Length;  PBFs.Initial();  PBFs.maptab=Maps;
                    memset(PBFs.labels,0,Length*sizeof(int));PBFs.Post_Type=atoi(SubSubType);break;
            case 'R':RBFs.k=ii;  RBFs.length=Length;  RBFs.Initial();  RBFs.maptab=Maps;  break;
            case 'T':TBFs.k=ii;  TBFs.length=Length;  TBFs.Initial();  TBFs.maptab=Maps;
                    TBFs.Thres_Type=atoi(SubSubType);  TBFs.In_Enumerate();  break;
            case 'U':UBFs.k=ii;  UBFs.length=Length;  UBFs.Initial();  UBFs.maptab=Maps;
                    UBFs.Unate_Type=atoi(SubSubType);  break;
            default: break;}
        for(jj=0;jj<Loop;jj++){
            switch(BoolFunType[0]){
                case 'C':CBFs.Cana_Gen();break; 
                case 'D':DBFs.Domin_Gen();break;
                case 'E':EBFs.Effec_Gen();break;
                case 'P':PBFs.Post_Gen(); break;
                case 'T':TBFs.Thres_Gen(0.5);break;
                case 'U':UBFs.Unate_Gen();break;
                case 'R':default:RBFs.RBF_p(0.5);break;}
            sensty_s[jj]=BFs.Sensitivity();// Sensitivity of normal SBFs 
            energy_s[jj]=BFs.Energy();// Energy of ...
            BFs.Exchange_Map();
            sensty_r[jj]=BFs.Sensitivity();// Sensitivity of exchanged SBFs 
            energy_r[jj]=BFs.Energy();// Energy of ...
            io[1]=BFs.Num_1();io[0]=Length-io[1];
            if(io[1]>io[0]){
                MajorBia[jj]=io[1]/Length;}
            else {
                MajorBia[jj]=io[0]/Length;}
            if(fabs(m_sen-sensty_s[jj])>1e-6){
                sensty__[jj]=(sensty_r[jj]-sensty_s[jj])/(m_sen-sensty_s[jj]);}
            else sensty__[jj]=0;
            if(fabs(energy_s[jj])>1e-6){
                energy__[jj]=(energy_s[jj]-fabs(energy_r[jj]))/energy_s[jj];}
            else energy__[jj]=0;}
        // Calculate Av and SD of each measure and output.
        fprintf(fpss,"%d,",ii);
        Mean_SD(sensty_s,Loop,io);fprintf(fpss,"%f,%f," ,io[0],io[1]);
        Mean_SD(sensty_r,Loop,io);fprintf(fpss,"%f,%f," ,io[0],io[1]);
        Mean_SD(sensty__,Loop,io);fprintf(fpss,"%f,%f," ,io[0],io[1]);
        Mean_SD(energy_s,Loop,io);fprintf(fpss,"%f,%f," ,io[0],io[1]);
        Mean_SD(energy_r,Loop,io);fprintf(fpss,"%f,%f," ,io[0],io[1]);
        Mean_SD(energy__,Loop,io);fprintf(fpss,"%f,%f," ,io[0],io[1]);
        Mean_SD(MajorBia,Loop,io);fprintf(fpss,"%f,%f\n",io[0],io[1]);
        // Clear and reset class-slot for next group of parameter.
        switch(BoolFunType[0]){
            case 'C':CBFs.Cana_Reset();break;
            case 'D':DBFs.Reset();break;case 'E':EBFs.Reset();break;
            case 'P':PBFs.Reset();break;case 'T':TBFs.Reset();break;
            case 'U':UBFs.Reset();break;case 'R':RBFs.Reset();break;}
        BFs.Reset();
        Length<<=1;
        free(Maps);}
    B_Sensitivity_Analysis::Delete();
    fclose(fpss);
}
// Sensitivity & Energy execute codes.
void exe_B_Sensitivity_Analysis(char **argv){
    B_Sensitivity_Analysis sensitivity;
    char type=argv[1][0];
    switch(type){
        case 'R':sensitivity.Subtype(argv[1],argv[2],argv[3],NULL,NULL);break;
        case 'C':sensitivity.Subtype(argv[1],argv[2],argv[3],argv[4],argv[5]);break;
        case 'P':sensitivity.Subtype(argv[1],argv[2],argv[3],argv[4],NULL);break;
        case 'D':sensitivity.Subtype(argv[1],argv[2],argv[3],argv[4],NULL);break;
        case 'E':sensitivity.Subtype(argv[1],argv[2],argv[3],NULL, NULL);break;
        case 'U':sensitivity.Subtype(argv[1],argv[2],argv[3],argv[4],NULL);break;
        case 'T':sensitivity.Subtype(argv[1],argv[2],argv[3],argv[4],NULL);break;
        default: printf("Error Type Arguments!\n");break;}
}
#endif // BOOLEAN_FUNCTION_SENSITIVITY_H_INCLUDED
