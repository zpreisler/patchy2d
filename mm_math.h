#ifndef MM_MATH_H
#define MM_MATH_H
#define NDEXP 32
#define SQR(x) ((x)*(x))
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
extern void boundary(__m128d *r,__m128d box);
extern __m128d _mm_dist(__m128d p,__m128d q,__m128d b);
extern __m128d _mm_dist_uy(__m128d p,__m128d q,__m128d b,double uy);
extern double length2(__m128d x);
extern double length(__m128d x);
extern __m128d rnd11();
extern __m128d sincosa(double a);
extern double dot(__m128d a,__m128d b);
extern __m128d normalize(__m128d a);
extern __m128d rot22(__m128d a,__m128d b01);
extern __m128d rot2w(__m128d a,double w);
extern unsigned int rndn(unsigned int n);
#endif
