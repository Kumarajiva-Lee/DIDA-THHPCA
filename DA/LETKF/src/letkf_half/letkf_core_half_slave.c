#include <stdio.h>
#include <slave.h>
#include <simd.h>
#include <crts.h>
#include <math.h>
#include "letkf_core_half.h"
#include "../../../utils/da_config.inc"

//scaling factor
#define alpha 1.0

void hxens_sqrt(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nobsl = info.nobsl; int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = 64;

    int bsize = (nobsl + threads - 1) / threads;
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nobsl);
    int length = letkf_min(bend - bstart, MAX_SPM / nanals / 2);

    yb_ctype hxens_sum[nanals];
    for (int i = 0; i < nanals; i++)
        hxens_sum[i] = 0;
    //if (id == 0) printf("%d %d %d\n", nanals, length, sizeof(_Float16));

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype hxens_in[length * nanals];
        yb_ctype hxens[length * nanals];
        float rdiaginv[length];
        yb_ctype rloc[length];
        ayb_ctype dep[length];

        CRTS_dma_get(hxens_in, info.hxens_in + bstart * nanals, sizeof(yb_ctype) * length * nanals);
        CRTS_dma_get(rloc, info.rloc + bstart, sizeof(yb_ctype) * length);
        CRTS_dma_get(rdiaginv, info.rdiaginv + bstart, sizeof(float) * length);
        CRTS_dma_get(dep, info.dep + bstart, sizeof(ayb_ctype) * length);
        //CRTS_dma_get(rloc, info.rloc + bstart, sizeof(yb_ctype) * length * nanals);

        for (int i = 0; i < length; i++)
        {
            rloc[i] = rloc[i] * rdiaginv[i];
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                hxens_sum[j] = hxens_sum[j] + rloc[i] * hxens_in[k] * dep[i];
            }
            //rloc[i] = sqrt(rloc[i]);
            rloc[i] = sqrt(rloc[i]);
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                //scaling
                hxens[k] = hxens_in[k] * rloc[i] / alpha;
                //double pnum = hxens[k];
                //if (id == 0) printf("%.5f %.5f %.5f\n", hxens_in[k], rloc[i], pnum);
            }
        }
        //if (id == 0) printf("\n");
        CRTS_dma_put(info.rrloc + bstart, rloc, sizeof(yb_ctype) * length);
        CRTS_dma_put((yb_ctype*)info.hxens + bstart * nanals, hxens, sizeof(yb_ctype) * length * nanals);
    }

    CRTS_dma_put(info.hxens_sum_64 + id * nanals, hxens_sum, sizeof(yb_ctype) * nanals);

    CRTS_ssync_array();

    bsize = (nanals + threads - 1) / threads;
    bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    length = letkf_min(bend - bstart, MAX_SPM / threads);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype hxens_sum_64[length * threads];
        yb_ctype hxens_sum[length];

        CRTS_dma_get_stride(hxens_sum_64, info.hxens_sum_64 + bstart, sizeof(yb_ctype) * length * threads, \
        sizeof(yb_ctype) * length, sizeof(yb_ctype) * (nanals - length));

        for (int i = 0; i < length; i++)
        {
            hxens_sum[i] = 0;
        }
        for (int i = 0; i < threads; i++)
        {
            for (int j = 0; j < length; j++)
            {
                int k = i * length + j;
                hxens_sum[j] += hxens_sum_64[k];
                //if (id == 0) printf("%8.4f ", ow_64[k]);
            }
        }
        //if (id == 0) printf("\n");

        CRTS_dma_put(info.hxens_sum_64 + bstart, hxens_sum, sizeof(yb_ctype) * length);
    }

}

void work1_half2double(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = THREAD_NUM;

    int bsize = (nanals + threads - 1) / threads;
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    int length = letkf_min(bend - bstart, MAX_SPM / nanals / 2);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype work1[length * nanals];
        //CRTS_dma_get(work1, info.work1 + bstart * nanals, sizeof(yb_ctype) * length * nanals);
        CRTS_dma_get(work1, (yb_ctype*)info.work1_half + bstart * nanals, sizeof(yb_ctype) * length * nanals);

        for (int i = 0; i < length; i++)
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                //work1[k] = work1_half[k];
                //scaling
                work1[k] = work1[k] * alpha * alpha;
            }
        
        CRTS_dma_put(info.work1 + bstart * nanals, work1, sizeof(yb_ctype) * length * nanals);
    }
}

void work1_double2half(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = THREAD_NUM;

    int bsize = (nanals + threads - 1) / threads;
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    int length = letkf_min(bend - bstart, MAX_SPM / nanals / 3);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype eivec[length * nanals];
        yb_ctype work1[length * nanals];
        yb_ctype eival[length];
        CRTS_dma_get(work1, info.work1 + bstart * nanals, sizeof(yb_ctype) * length * nanals);
        CRTS_dma_get(eival, info.eival + bstart, sizeof(yb_ctype) * length);

        for (int i = 0; i < length; i++)
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                eivec[k] = work1[k];
                work1[k] = work1[k] / eival[i];
            }

        CRTS_dma_put(info.work1 + bstart * nanals, work1, sizeof(yb_ctype) * length * nanals);
        CRTS_dma_put(info.eivec_double + bstart * nanals, eivec, sizeof(yb_ctype) * length * nanals);
 
    }
}

void times_rrloc(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nobsl = info.nobsl; int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = THREAD_NUM;

    int bsize = (nobsl + threads - 1) / threads;
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nobsl);
    int length = letkf_min(bend - bstart, MAX_SPM / nanals / 2);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype hxens[length * nanals];
        yb_ctype rrloc[length];
        CRTS_dma_get(hxens, info.hxens + bstart * nanals, sizeof(yb_ctype) * length * nanals);
        CRTS_dma_get(rrloc, info.rrloc + bstart, sizeof(yb_ctype) * length);

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                hxens[k] = hxens[k] * rrloc[i];
            }
        }

        CRTS_dma_put(info.hxens + bstart * nanals, hxens, sizeof(yb_ctype) * length * nanals);

    }

}

void ow_calc(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nobsl = info.nobsl; int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = THREAD_NUM;

    int bsize = (nanals + threads - 1) / threads;
    //bsize = letkf_max(bsize, 2);
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    int length = letkf_min(bend - bstart, MAX_SPM / nanals);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype pa_double[length * nanals];
        yb_ctype hxens_sum_64[nanals];
        yb_ctype ow[length];

        //CRTS_dma_get_stride(pa, (_Float16*)info.pa + bstart, sizeof(_Float16) * length * nanals, \
        sizeof(_Float16*) * length, sizeof(_Float16*) * (nanals - length));
        CRTS_dma_get(pa_double, info.pa_double + bstart * nanals, sizeof(yb_ctype) * length  * nanals);
        CRTS_dma_get(hxens_sum_64, info.hxens_sum_64, sizeof(yb_ctype) * nanals);

        for (int i = 0; i < length; i++)
        {
            ow[i] = 0;
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                ow[i] = ow[i] + pa_double[k] * hxens_sum_64[j];
                //if (id == 31) printf("%16.12f %8.4f\n", (double)pa_double[k], hxens_sum_64[j]);
            }
        }
        // if (id == 31)
        // {
        //     printf("hxens_sum_64\n");
        //     for (int i = 0; i < nanals; i++)
        //     {
        //         printf("%8.4f ", hxens_sum_64[i]);
        //     }
        //     printf("\n");
        //     printf("ow %8.4f\n", ow[0]);
        // }

        CRTS_dma_put(info.overline_w + bstart, ow, sizeof(yb_ctype) * length);
    }

    // CRTS_ssync_array();

    // bsize = (nanals + threads - 1) / threads;
    // bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    // length = letkf_min(bend - bstart, MAX_SPM / threads);

    // for (bstart; bstart < bend; bstart += length)
    // {
    //     length = letkf_min(bend - bstart, length);

    //     yb_ctype ow_64[length * threads];
    //     yb_ctype ow[length];

    //     CRTS_dma_get_stride(ow_64, info.ow_64 + bstart, sizeof(yb_ctype) * length * threads, \
    //     sizeof(yb_ctype) * length, sizeof(yb_ctype) * (nanals - length));

    //     for (int i = 0; i < length; i++)
    //     {
    //         ow[i] = 0;
    //     }
    //     for (int i = 0; i < threads; i++)
    //     {
    //         for (int j = 0; j < length; j++)
    //         {
    //             int k = i * length + j;
    //             ow[j] += ow_64[k];
    //             //if (id == 0) printf("%8.4f ", ow_64[k]);
    //         }
    //     }
    //     //if (id == 0) printf("\n");

    //     CRTS_dma_put(info.overline_w + bstart, ow, sizeof(yb_ctype) * length);
    // }

}

void sqrt_pa(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = THREAD_NUM;

    int bsize = (nanals + threads - 1) / threads;
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    int length = letkf_min(bend - bstart, MAX_SPM / 3 / nanals);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        yb_ctype eivec_double[nanals * length];
        yb_ctype work1[nanals * length];
        yb_ctype eival[length];

        CRTS_dma_get(eivec_double, info.eivec_double + bstart * nanals, sizeof(yb_ctype) * length * nanals);
        CRTS_dma_get(eival, info.eival + bstart, sizeof(yb_ctype) * length);

        for (int i = 0; i < length; i++)
        {
            yb_ctype rho = sqrt( (nanals - 1) / eival[i] );
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                //if (id == 0) printf("%8.4f ", (double)work1_half[k]);
                work1[k] = eivec_double[k] * rho;
            }
            //if (id == 0) printf("end id1\n");
        }

        CRTS_dma_put((yb_ctype*)info.work1 + bstart * nanals, work1, sizeof(yb_ctype) * length * nanals);
        //CRTS_dma_put((yb_ctype)info.eivec + bstart * nanals, eivec, sizeof(_Float16) * length * nanals);
    }

}

void trans_half2double(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals;

    int id = CRTS_smng_get_tid();
    int threads = THREAD_NUM;

    int bsize = (nanals + threads - 1) / threads;
    int bstart = bsize * id, bend = letkf_min(bsize * (id + 1), nanals);
    int length = letkf_min(bend - bstart, MAX_SPM / nanals / 2);

    yb_ctype ow[nanals];
    CRTS_dma_get(ow, info.overline_w, sizeof(yb_ctype) * nanals);

    for (bstart; bstart < bend; bstart += length)
    {
        length = letkf_min(bend - bstart, length);

        //_Float16 trans_half[nanals * length];
        yb_ctype trans[nanals * length];

        CRTS_dma_get(trans, (yb_ctype*)info.trans + bstart * nanals, sizeof(yb_ctype) * length * nanals);

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                trans[k] = trans[k] + ow[j];
            }
        }

        CRTS_dma_put(info.trans + bstart * nanals, trans, sizeof(yb_ctype) * length * nanals);
    }
}

void debug_show(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals; int nobsl = info.nobsl;

    if (CRTS_smng_get_tid() == 0)
    {   
        _Float16* pa = info.trans_half;
        for (int i = 0; i < nanals; i++)
        {
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                double tmp = pa[k];
                printf("%8.4f ", tmp);
            }
            printf("\n");
        }
    }
}

void debug_show2(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals; int nobsl = info.nobsl;

    if (CRTS_smng_get_tid() == 0)
    {   
        _Float16* pa = info.eivec;
        for (int i = 0; i < nanals; i++)
        {
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                double tmp = pa[k];
                printf("%8.4f ", tmp);
            }
            printf("\n");
        }
    }
}

void debug_show3(void* _info)
{
    half_info info;
    CRTS_dma_get(&info, _info, sizeof(half_info));
    int nanals = info.nanals; int nobsl = info.nobsl;

    if (CRTS_smng_get_tid() == 0)
    {   
        _Float16* pa = info.work1_half;
        for (int i = 0; i < nanals; i++)
        {
            for (int j = 0; j < nanals; j++)
            {
                int k = i * nanals + j;
                double tmp = pa[k];
                printf("%8.4f ", tmp);
            }
            printf("\n");
        }
    }
}
