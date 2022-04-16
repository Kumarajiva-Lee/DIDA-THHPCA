#define CYCLIC_TRIDIAGONAL_HOST_SIMD
//#define CYCLIC_TRIDIAGONAL_SLAVE_SIMD
#ifdef CYCLIC_TRIDIAGONAL_TEST
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#endif
#if defined(CYCLIC_TRIDIAGONAL_HOST_SIMD) || defined(CYCLIC_TRIDIAGONAL_SLAVE_SIMD)
#include<simd.h>
#endif

#define RSTR __restrict__

/* loop unrolling */
template<typename T, typename Func>
inline void for_each(T (&x)[1], Func && func) {
    func(x[0]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[2], Func && func) {
    func(x[0]);
    func(x[1]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[8], Func && func) {
    func(x[0]);
    func(x[1]);
    func(x[2]);
    func(x[3]);
    func(x[4]);
    func(x[5]);
    func(x[6]);
    func(x[7]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[1], T (&y)[1], Func && func) {
    func(x[0], y[0]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[2], T (&y)[2], Func && func) {
    func(x[0], y[0]);
    func(x[1], y[1]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[8], T (&y)[8], Func && func) {
    func(x[0], y[0]);
    func(x[1], y[1]);
    func(x[2], y[2]);
    func(x[3], y[3]);
    func(x[4], y[4]);
    func(x[5], y[5]);
    func(x[6], y[6]);
    func(x[7], y[7]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[1], T (&y)[1], T (&z)[1], Func && func) {
    func(x[0], y[0], z[0]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[2], T (&y)[2], T (&z)[2], Func && func) {
    func(x[0], y[0], z[0]);
    func(x[1], y[1], z[1]);
}

template<typename T, typename Func>
inline void for_each(T (&x)[8], T (&y)[8], T (&z)[8], Func && func) {
    func(x[0], y[0], z[0]);
    func(x[1], y[1], z[1]);
    func(x[2], y[2], z[2]);
    func(x[3], y[3], z[3]);
    func(x[4], y[4], z[4]);
    func(x[5], y[5], z[5]);
    func(x[6], y[6], z[6]);
    func(x[7], y[7], z[7]);
}


/*
 * on input:
 *     a is the lower sub-diagonal.
 *     b is the main diagonal.
 *     c is the upper sub-diagonal.
 *     rhs is the right-hand-side vector.
 *     a[0] is the upper-right element.
 *     c[n - 1] is the lower-left element.
 * 
 * on output:
 *     a and c are overwritten by internal data. That means, a and c are destroyed
 *     rhs will be overwritten by the solution.
 */
template<typename val_t, typename RHS, size_t nvec>
void cyclic_tridiagonal(int n, val_t* RSTR a, const val_t* RSTR b, val_t* RSTR c, RHS (*RSTR rhs)[nvec]) {
    int i;
    val_t un2, bi, v = c[n - 1], r = b[n - 1];
    RHS d[nvec], cons1, cons2;
    // d = rhs[n - 1];
    for_each(d, rhs[n - 1], [](RHS& x, RHS y) { x = y; } );
    a[0]   = a[0]   / b[0];
    c[0]   = c[0]   / b[0];
    // rhs[0] = rhs[0] / b[0];
    cons1 = b[0];
    for_each(rhs[0], [=](RHS& x) { x = x / cons1; } );
    for(i = 1; i < n - 1; ++ i) {
        // row [i - 1] eliminate row[i]
        bi     = (b[i] - c[i - 1] * a[i]);
        c[i]   = c[i] / bi;

        // rhs[i] = (rhs[i] - rhs[i - 1] * a[i]) / bi;
        cons1 = a[i];
        cons2 = bi;
        for_each(rhs[i], rhs[i - 1], [=](RHS& x, RHS y) { x = (x - y * cons1) / cons2; } );
        a[i]   = ( - a[i - 1] * a[i]) / bi;  // a[i - 1] has been overwritten by u[i - 1]
        // row [i - 1] eliminate row[n - 1]
        // d = d - rhs[i - 1] * v;
        cons1 = v;
        for_each(d, rhs[i - 1], [=](RHS& x, RHS y) { x = x - y * cons1; } );
        r = r -   a[i - 1] * v;
        v =   -   c[i - 1] * v;
    }
    // alias
    v = v + a[n - 1];
    un2 = a[n - 2] + c[n - 2];
    // solve rhs[n - 1]
    bi = r - v * un2;
    // rhs[n - 1] = (d - v * rhs[n - 2]) / bi;
    cons1 = v;
    cons2 = bi;
    for_each(rhs[n - 1], d, rhs[n - 2], [=](RHS& x, RHS y, RHS z) { x = (y - cons1 * z) / cons2; } );

    // solve rhs[n - 2]
    // rhs[n - 2] = rhs[n - 2] - un2 * rhs[n - 1];
    cons1 = un2;
    for_each(rhs[n - 2], rhs[n - 1], [=](RHS& x, RHS y) { x = x - cons1 * y; });

    // solve rhs
    for(i = n - 3; i >= 0; --i) {
        // rhs[i] = rhs[i] - c[i] * rhs[i + 1] - a[i] * rhs[n - 1];
        cons1 = c[i];
        cons2 = a[i];
        for_each(rhs[i], rhs[i + 1], rhs[n - 1], [=](RHS& x, RHS y, RHS z) { x = x - cons1 * y - cons2 * z; } );
    }
}

// Fortran callable
extern "C" void cyclic_tridiagonal_d_(const int* n, double* RSTR a, const double* RSTR b, double* RSTR c, double* RSTR rhs) {
    typedef double (*ptr)[1];
    cyclic_tridiagonal(*n, a, b, c, reinterpret_cast<ptr>(rhs));
}


extern "C" void cyclic_tridiagonal_d8_(const int* n, double* RSTR a, const double* RSTR b, double* RSTR c, double (*RSTR rhs)[8]) {
#if defined(CYCLIC_TRIDIAGONAL_HOST_SIMD)
    typedef doublev4 (*ptr)[2];
#elif defined(CYCLIC_TRIDIAGONAL_SLAVE_SIMD)
    typedef doublev8 (*ptr)[1];
#else
    typedef double (*ptr)[8];
#endif
    cyclic_tridiagonal(*n, a, b, c, reinterpret_cast<ptr>(rhs));
}

#ifdef CYCLIC_TRIDIAGONAL_TEST
#define SIZE 2400

int main() {
    double em, ea;
    double a[SIZE], b[SIZE], c[SIZE], x[SIZE], rhs[SIZE];
    alignas(64) double xv[SIZE][8];
    alignas(64) double rhsv[SIZE][8];
    int n = SIZE;
    for(int i = 0; i < n; ++i) {
        a[i] = (double)rand() / (double)(RAND_MAX / 2) - 1.0;
        b[i] = (double)rand() / (double)(RAND_MAX / 2) - 1.0; // better be diagonal dominant
        c[i] = (double)rand() / (double)(RAND_MAX / 2) - 1.0;
        x[i] = (double)rand() / (double)(RAND_MAX / 2) - 1.0;
    }
    rhs[0] = a[0] * x[n-1] + b[0] * x[0] + c[0] * x[1];
    for(int i = 1; i < n - 1; ++i) {
        rhs[i] = a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1];
    }
    rhs[n-1] = a[n-1] * x[n-2] + b[n-1] * x[n-1] + c[n-1] * x[0];
    cyclic_tridiagonal_d_(&n, a, b, c, rhs);

    em = 0.0;
    ea = 0.0;
    for(int i = 0; i < n; ++i) {
        double e = fabs(x[i] - rhs[i]);
        em = (em < e) ? e : em;
        ea += e * e;
    }
    printf("1 right-hand-side vector:\n\tmaximum error = %e\n\taverage error = %e\n", em, sqrt(ea / n));
    
    for(int i = 0; i < n; ++i) {
        a[i] = (double)rand() / (double)(RAND_MAX / 2) - 1.0;
        c[i] = (double)rand() / (double)(RAND_MAX / 2) - 1.0;
        for(int j = 0; j < 8; ++j) {
            xv[i][j] = (double)rand() / (double)(RAND_MAX / 2) - 1.0;
        }
    }
    for(int j = 0; j < 8; ++j) {
        rhsv[0][j] = a[0] * xv[n-1][j] + b[0] * xv[0][j] + c[0] * xv[1][j];
    }
    for(int i = 1; i < n - 1; ++i) {
        for(int j = 0; j < 8; ++j) {
            rhsv[i][j] = a[i] * xv[i-1][j] + b[i] * xv[i][j] + c[i] * xv[i+1][j];
        }
    }
    for(int j = 0; j < 8; ++j) {
        rhsv[n-1][j] = a[n-1] * xv[n-2][j] + b[n-1] * xv[n-1][j] + c[n-1] * xv[0][j];
    }

    cyclic_tridiagonal_d8_(&n, a, b, c, rhsv);
    
    em = 0.0;
    ea = 0.0;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < 8; ++j) {
            double e = fabs(xv[i][j] - rhsv[i][j]);
            em = (em < e) ? e : em;
            ea += e * e;
        }
    }
    printf("8 right-hand-side vectors:\n\tmaximum error = %e\n\taverage error = %e\n", em, sqrt(ea / n));
}

#endif
