#ifndef LETKF_CORE_HALF_
#define LETKF_CORE_HALF_

#include "../../../utils/da_config.inc"

#if (yb_size == 8)
    #define yb_ctype double
#else
    #define yb_ctype float
#endif

#if (ayb_size == 8)
    #define ayb_ctype double
#else
    #define ayb_ctype float
#endif

#define letkf_max(a, b) (((a)>(b))?(a):(b))
#define letkf_min(a, b) (((a)<(b))?(a):(b))

#define MAX_SPM (16*1024)
#define THREAD_NUM 64

typedef struct
{
    int nanals;
    int nobsl;
    const yb_ctype *hxens_in;
    void *hxens;
    const float *rdiaginv;
    const ayb_ctype *dep;
    const yb_ctype *rloc;
    yb_ctype *rrloc;
    yb_ctype *trans;
    void *trans_half;
    yb_ctype *overline_w;
    yb_ctype *ow_64;
    void *work1_half;
    yb_ctype *work1;
    void *eivec;
    yb_ctype *eivec_double;
    yb_ctype *eival;
    void *work2;
    void *pa;
    yb_ctype *pa_double;
    yb_ctype *hxens_sum_64;
} half_info;

#endif
