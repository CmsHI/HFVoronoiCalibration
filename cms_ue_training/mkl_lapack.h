#ifndef _MKL_LAPACK_H_
#define _MKL_LAPACK_H_

typedef int MKL_INT;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void dgelss_( const MKL_INT* m, const MKL_INT* n, const MKL_INT* nrhs, 
              double* a, const MKL_INT* lda, double* b, const MKL_INT* ldb,
              double* s, const double* rcond, MKL_INT* rank, double* work,
              const MKL_INT* lwork, MKL_INT* info );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MKL_LAPACK_H_ */
