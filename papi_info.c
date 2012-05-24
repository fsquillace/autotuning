
#include <papi.h>



int main(){
    int events[107];
events[0]=PAPI_VEC_DP;
events[1]=PAPI_L1_DCM;
events[2]=PAPI_L1_ICM;
events[3]=PAPI_L2_DCM;
events[4]=PAPI_L2_ICM;
events[5]=PAPI_L3_DCM;
events[6]=PAPI_L3_ICM;
events[7]=PAPI_L1_TCM;
events[8]=PAPI_L2_TCM;
events[9]=PAPI_L3_TCM;
events[10]=PAPI_CA_SNP;
events[11]=PAPI_CA_SHR;
events[12]=PAPI_CA_CLN;
events[13]=PAPI_CA_INV;
events[14]=PAPI_CA_ITV;
events[15]=PAPI_L3_LDM;
events[16]=PAPI_L3_STM;
events[17]=PAPI_BRU_IDL;
events[18]=PAPI_FXU_IDL;
events[19]=PAPI_FPU_IDL;
events[20]=PAPI_LSU_IDL;
events[21]=PAPI_TLB_DM;
events[22]=PAPI_TLB_IM;
events[23]=PAPI_TLB_TL;
events[24]=PAPI_L1_LDM;
events[25]=PAPI_L1_STM;
events[26]=PAPI_L2_LDM;
events[27]=PAPI_L2_STM;
events[28]=PAPI_BTAC_M;
events[29]=PAPI_PRF_DM;
events[30]=PAPI_L3_DCH;
events[31]=PAPI_TLB_SD;
events[32]=PAPI_CSR_FAL;
events[33]=PAPI_CSR_SUC;
events[34]=PAPI_CSR_TOT;
events[35]=PAPI_MEM_SCY;
events[36]=PAPI_MEM_RCY;
events[37]=PAPI_MEM_WCY;
events[38]=PAPI_STL_ICY;
events[39]=PAPI_FUL_ICY;
events[40]=PAPI_STL_CCY;
events[41]=PAPI_FUL_CCY;
events[42]=PAPI_HW_INT;
events[43]=PAPI_BR_UCN;
events[44]=PAPI_BR_CN;
events[45]=PAPI_BR_TKN;
events[46]=PAPI_BR_NTK;
events[47]=PAPI_BR_MSP;
events[48]=PAPI_BR_PRC;
events[49]=PAPI_FMA_INS;
events[50]=PAPI_TOT_IIS;
events[51]=PAPI_TOT_INS;
events[52]=PAPI_INT_INS;
events[53]=PAPI_FP_INS;
events[54]=PAPI_LD_INS;
events[55]=PAPI_SR_INS;
events[56]=PAPI_BR_INS;
events[57]=PAPI_VEC_INS;
events[58]=PAPI_RES_STL;
events[59]=PAPI_FP_STAL;
events[60]=PAPI_TOT_CYC;
events[61]=PAPI_LST_INS;
events[62]=PAPI_SYC_INS;
events[63]=PAPI_L1_DCH;
events[64]=PAPI_L2_DCH;
events[65]=PAPI_L1_DCA;
events[66]=PAPI_L2_DCA;
events[67]=PAPI_L3_DCA;
events[68]=PAPI_L1_DCR;
events[69]=PAPI_L2_DCR;
events[70]=PAPI_L3_DCR;
events[71]=PAPI_L1_DCW;
events[72]=PAPI_L2_DCW;
events[73]=PAPI_L3_DCW;
events[74]=PAPI_L1_ICH;
events[75]=PAPI_L2_ICH;
events[76]=PAPI_L3_ICH;
events[77]=PAPI_L1_ICA;
events[78]=PAPI_L2_ICA;
events[79]=PAPI_L3_ICA;
events[80]=PAPI_L1_ICR;
events[81]=PAPI_L2_ICR;
events[82]=PAPI_L3_ICR;
events[83]=PAPI_L1_ICW;
events[84]=PAPI_L2_ICW;
events[85]=PAPI_L3_ICW;
events[86]=PAPI_L1_TCH;
events[87]=PAPI_L2_TCH;
events[88]=PAPI_L3_TCH;
events[89]=PAPI_L1_TCA;
events[90]=PAPI_L2_TCA;
events[91]=PAPI_L3_TCA;
events[92]=PAPI_L1_TCR;
events[93]=PAPI_L2_TCR;
events[94]=PAPI_L3_TCR;
events[95]=PAPI_L1_TCW;
events[96]=PAPI_L2_TCW;
events[97]=PAPI_L3_TCW;
events[98]=PAPI_FML_INS;
events[99]=PAPI_FAD_INS;
events[100]=PAPI_FDV_INS;
events[101]=PAPI_FSQ_INS;
events[102]=PAPI_FNV_INS;
events[103]=PAPI_FP_OPS;
events[104]=PAPI_SP_OPS;
events[105]=PAPI_DP_OPS;
events[106]=PAPI_VEC_SP;


    int retval;
    // Initialize the library
    retval = PAPI_library_init(PAPI_VER_CURRENT);

    if (retval != PAPI_VER_CURRENT) {
       printf("PAPI library init error!\n");
          exit(1);
    }

    PAPI_event_info_t info;
    for(int i=0; i<107; i++){
        if (PAPI_query_event(events[i]) != PAPI_OK) {
            printf("No instruction counter? How lame.\n");
        }
        else{
            PAPI_get_event_info(events[i], &info);
            printf("************ %s **********\n", info.symbol);
            printf("symbol=%s\n", info.symbol);
            printf("short_descr=%s\n", info.short_descr);
            printf("long_descr=%s\n", info.long_descr);
            printf("derived=%s\n", info.derived);
            printf("postfix=%s\n", info.postfix);
            printf("code=%d\n", info.code);
            printf("name=%s\n", info.name);
            printf("note=%s\n\n", info.note);
        }

    }


    return 0;
}
