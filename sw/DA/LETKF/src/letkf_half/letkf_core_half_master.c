//#include "../../utils/da_config.inc"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <crts.h>
#include <swhgemm.h>
#include <gptl.h>
#include "letkf_core_half.h"
#include "../../../utils/da_config.inc"
#define DEBUG 1

extern void SLAVE_FUN(hxens_sqrt)(void *);
extern void SLAVE_FUN(work1_half2double)(void *);
extern void SLAVE_FUN(work1_double2half)(void *);
extern void SLAVE_FUN(times_rrloc)(void *);
extern void SLAVE_FUN(ow_calc)(void *);
extern void SLAVE_FUN(sqrt_pa)(void *);
extern void SLAVE_FUN(trans_half2double)(void *);
extern void SLAVE_FUN(debug_show)(void *);
extern void SLAVE_FUN(debug_show2)(void *);
extern void SLAVE_FUN(debug_show3)(void *);

void letkf_core_half_c_(const int *nanals_in, const int *nobsl_in, const yb_ctype *hxens_in, const float *rdiaginv, \
const ayb_ctype *dep, const yb_ctype *rloc, yb_ctype *trans, yb_ctype *overline_w, \
const double *inflation_factor, const int *pflag)
{

    int nanals = *nanals_in, nobsl = *nobsl_in;
    int nob, nanal;
    int nobsl2 = nobsl;
    if (nobsl2 % 2 == 1) nobsl2 += 1;

    __attribute__ ((aligned(64))) yb_ctype hxens[nobsl2 * nanals];
    __attribute__ ((aligned(64))) yb_ctype rrloc[nobsl2];

    __attribute__ ((aligned(64))) int work1_half[nanals * nanals];
    __attribute__ ((aligned(64))) yb_ctype work1[nanals * nanals];

    __attribute__ ((aligned(64))) int eivec[nanals * nanals];
    __attribute__ ((aligned(64))) yb_ctype eivec_double[nanals * nanals];
    __attribute__ ((aligned(64))) yb_ctype eival[nanals];

    __attribute__ ((aligned(64))) int pa[nanals * nanals];
    __attribute__ ((aligned(64))) yb_ctype pa_double[nanals * nanals];

    __attribute__ ((aligned(64))) int work2[nobsl2 * nanals];

    __attribute__ ((aligned(64))) yb_ctype ow_64[THREAD_NUM * nanals];

    __attribute__ ((aligned(64))) int trans_half[nanals * nanals];

    __attribute__ ((aligned(64))) yb_ctype hxens_sum_64[THREAD_NUM * nanals];

    int sizework = 2*nanals*nanals+8*nanals, sizeiwork = 6*nanals;
    __attribute__ ((aligned(64))) yb_ctype work[sizework];
    __attribute__ ((aligned(64))) int iwork[sizeiwork];

    half_info info;
    info.nanals = nanals; info.nobsl = nobsl;
    info.hxens = hxens; info.hxens_in = hxens_in;
    info.rdiaginv = rdiaginv; info.dep = dep;
    info.rloc = rloc; info.rrloc = rrloc;
    info.trans = trans; info.trans_half = trans_half;
    info.overline_w = overline_w; info.ow_64 = ow_64;
    info.work1 = work1; info.work1_half = work1_half;
    info.eivec_double = eivec_double; info.eivec = eivec;
    info.eival = eival;
    info.work2 = work2;
    info.pa = pa; info.pa_double = pa_double;
    info.hxens_sum_64 = hxens_sum_64;

    for (int i = 0; i < THREAD_NUM * nanals; i++) ow_64[i] = 0;

    CRTS_init();

    CRTS_athread_spawn(hxens_sqrt, &info);
    CRTS_athread_join();

    // for (int i = 0; i < nanals; i++)
    // {
    //     printf("%8.4f ", hxens_sum_64[i]);
    // }
    // printf("\n");

    float alpha = 1.0, beta = 0.0;
    double alpha2 = 1.0, beta2 = 0.0;
    // SW_TRANSPOSE TA = swNoTrans, TB = swTrans;
    // hgemm_(&TA, &TB, &nanals, &nanals, &nobsl, &alpha, hxens, &nanals, \
    //         hxens, &nobsl, &beta, work1_half, &nanals);
#if (DEBUG == 1)
    GPTLstart("big gemm");
#endif
    char TA = 'N', TB = 'T';
    dgemm_(&TA, &TB, &nanals, &nanals, &nobsl, &alpha2, hxens, &nanals, \
            hxens, &nanals, &beta2, work1, &nanals);
#if (DEBUG == 1)
    GPTLstop("big gemm");
#endif
    
    //printf("hgemm done!\n");

    // CRTS_athread_spawn(work1_half2double, &info);
    // CRTS_athread_join();

    //for (int i = 0; i < nanals; i++)
    //    work1[i * nanals + i] += (nanals - 1) / (*inflation_factor);

	if (*pflag)
	{
		printf("work1\n");
		for (int i = 0; i < nanals; i++)
		{
			for (int j = 0; j < nanals; j++)
			{
     	        int k = i * nanals + j;
     	        printf("%.4f ", work1[k]);
     	    }
			printf("\n");
		}
	}

    // if (*pflag)
	// {
    //     printf("work1 for svd\n");
	// 	for (int i = 0; i < nanals; i++)
	// 	{
	// 		for (int j = 0; j < nanals; j++)
	// 		{
    //  	        int k = i * nanals + j;
    //  	        printf("%.4f ", work1[k]);
    //  	    }
	// 		printf("\n");
	// 	}
    // }

    char jobz = 'V', uplo = 'L';
    int ierr;

#if (DEBUG == 1)
   GPTLstart("eig decomposition");
#endif
    dsyevd_(&jobz, &uplo, &nanals, work1, &nanals, eival, work, &sizework, iwork, &sizeiwork, &ierr);
#if (DEBUG == 1)
   GPTLstop("eig decomposition");
#endif
    //CRTS_athread_join();
	
	// if (*pflag)
	// {
    //     printf("eival ori!\n");
    //     for (int i = 0; i < nanals; i++)
    //         printf("%.4f ", eival[i]);
	// 	printf("\n");
	// }

	if (eival[0] < 0)
	{
		yb_ctype eival_min = eival[0];
		for (int i = 0; i < nanals; i++)
			eival[i] = eival[i] - eival_min;
	}
	for (int i = 0; i < nanals; i++)
		eival[i] = eival[i] + (nanals - 1) / (*inflation_factor);

	if (*pflag)
	{
        printf("eival!\n");
        for (int i = 0; i < nanals; i++)
            printf("%.4f ", eival[i]);
        printf("\n");
        printf("eivec!\n");
        for (int i = 0; i < nanals * nanals; i++)
			printf("%.4f ", work1[i]);
		printf("\n");
	}

    //if (ierr != 0 && *pflag)
    if (ierr != 0)
    {
        // printf("error nobsl = %d\n", nobsl);
        // printf("hxens_sum_64!\n");
        // for (int i = 0; i < nanals; i++)
        //     printf("%.4f ", hxens_sum_64[i]);
        // printf("\n");
        // printf("work1_half\n");
        // CRTS_athread_spawn(debug_show2, &info);
        // CRTS_athread_join();
        // printf("hxens\n");
        // CRTS_athread_spawn(debug_show, &info);
        // CRTS_athread_join();
        // for (int i = 0; i < nanals; i++)
		// {
		// 	// for (int j = 0; j < nanals; j++)
		// 	// {
     	//     //     int k = i * nanals + j;
     	//     //     printf("%.4f ", hxens_in[k]);
     	//     // }
		// 	// printf("\n");
		// }
        printf("dsyev failed, ierr=%d\n", ierr);
        //mpi_abort();
    }

    CRTS_athread_spawn(work1_double2half, &info);
    CRTS_athread_join();

#if (DEBUG == 1)
   GPTLstart("gemm");
#endif
    TA = 'N'; TB = 'T';
    //hgemm_(&TA, &TB, &nanals, &nanals, &nanals, &alpha, work1_half, &nanals, \
            eivec, &nanals, &beta, pa, &nanals);
    dgemm_(&TA, &TB, &nanals, &nanals, &nanals, &alpha2, work1, &nanals, \
            eivec_double, &nanals, &beta2, pa_double, &nanals);
#if (DEBUG == 1)
    GPTLstop("gemm");
#endif
	
	// if (*pflag)
	// {
	// 	printf("pa_double\n");
	// 	for (int i = 0; i < nanals; i++)
	// 	{
	// 		for (int j = 0; j < nanals; j++)
	// 		{
    //  	        int k = i * nanals + j;
    //  	        printf("%.5f ", pa_double[k]);
    //  	    }
	// 		printf("\n");
	// 	}
	// }

	// if (*pflag)
	// {
	// 	printf("pa\n");
	// 	CRTS_athread_spawn(debug_show2, &info);
	//     CRTS_athread_join();
	// }
	// if (*pflag)
	// {
	// 	printf("pa\n");
	// 	for (int i = 0; i < nanals; i++)
	// 	{
	// 		for (int j = 0; j < nanals; j++)
	// 		{
    //  	        int k = i * nanals + j;
    //  	        printf("%.5f ", pa[k]);
    //  	    }
	// 		printf("\n");
	// 	}
	// }
    
    //printf("hgemm done!\n");

    // CRTS_athread_spawn(times_rrloc, &info);
    // CRTS_athread_join();

    // //printf("times_rrloc done!\n");

    // TA = 'N'; TB = 'N';
    // hgemm_(&TA, &TB, &nanals, &nobsl, &nanals, &alpha, pa, &nanals, \
    //         hxens, &nanals, &beta, work2, &nanals);

    //printf("hgemm done!\n");

    CRTS_athread_spawn(ow_calc, &info);
    CRTS_athread_join();

    // printf("ow_64\n");
    // for (int i = 0; i < THREAD_NUM; i++)
    // {
    //     for (int j = 0; j < nanals; j++)
    //     {
    //         int k = i * nanals + j;
    //         printf("%8.4f ", ow_64[k]);
    //     }
    //     printf("\n");
    // }

    // printf("ow1 ");
    // for (int i = 0; i < nanals; i++)
    // {
    //     double tmp = 0;
    //     for (int j = 0; j < THREAD_NUM; j++)
    //     {
    //         int k = j * nanals + i;
    //         tmp = tmp + ow_64[k];
    //         //printf("%.4f ", ow_64[k]);
    //     }
    //     printf("%.4f ", tmp);
    // }
    // printf("\n");
	if (*pflag)
	{
		printf("ow ");
		for (int i = 0; i < nanals; i++)
    	{
    	    printf("%.4f ", overline_w[i]);
    	}
		printf("\n");
	}
    //printf("ow calc done!\n");

    CRTS_athread_spawn(sqrt_pa, &info);
    CRTS_athread_join();

    //printf("sqrt_pa done!\n");
    
#if (DEBUG == 1)
   GPTLstart("gemm");
#endif
    TA = 'N'; TB = 'T';
    dgemm_(&TA, &TB, &nanals, &nanals, &nanals, &alpha2, work1, &nanals, \
            eivec_double, &nanals, &beta2, trans, &nanals);
#if (DEBUG == 1)
   GPTLstop("gemm");
#endif
    
    // if (*pflag)
    // {
    //     printf("eivec_double\n");
    //     for (int i = 0; i < nanals; i++)
    //     {
    //         for (int j = 0; j < nanals; j++)
    //         {
    //             int k = i * nanals + j;
    //             printf("%8.4f ", eivec_double[k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("work1\n");
    //     for (int i = 0; i < nanals; i++)
    //     {
    //         for (int j = 0; j < nanals; j++)
    //         {
    //             int k = i * nanals + j;
    //             printf("%8.4f ", work1[k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("trans\n");
    //     for (int i = 0; i < nanals; i++)
    //     {
    //         for (int j = 0; j < nanals; j++)
    //         {
    //             int k = i * nanals + j;
    //             printf("%8.4f ", trans[k]);
    //         }
    //         printf("\n");
    //     }
    // }

	// if (*pflag)
	// {
    //     printf("trans_half\n");
	//     CRTS_athread_spawn(debug_show, &info);
	//     CRTS_athread_join();
	// }
	// if (*pflag)
	// {
    //     printf("eivec\n");
	//     CRTS_athread_spawn(debug_show2, &info);
	//     CRTS_athread_join();
	// }
	// if (*pflag)
	// {
    //     printf("eivec\n");
	//     CRTS_athread_spawn(debug_show3, &info);
	//     CRTS_athread_join();
	// }
    
    CRTS_athread_spawn(trans_half2double, &info);
    CRTS_athread_join();

    //printf("trans_half2double done!\n");

    if (*pflag)
	{
		printf("trans\n");
		for (int i = 0; i < nanals; i++)
		{
			for (int j = 0; j < nanals; j++)
			{
     	        int k = i * nanals + j;
     	        printf("%.5f ", trans[k]);
     	    }
			printf("\n");
		}
	}

    //void hgemm_(const char *TransA, const char *TransB, const int *M, const int *N, const int *K,
    //            const float *alpha, const void *A, const int *lda, const void *B, const int *ldb,
    //            const float beta, short *C, const int ldc);


}
