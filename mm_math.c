#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include <utils.h>
#include "dSFMT.h"
#include "mm_math.h"
extern dsfmt_t dsfmt;
void boundary(__m128d *r,__m128d box){
	__m128d t=_mm_div_pd(*r,box);	
	t=_mm_round_pd(t,0x1);
	t=_mm_mul_pd(t,box);
	*r=_mm_sub_pd(*r,t);
}
__m128d _mm_dist(__m128d p,__m128d q,__m128d b){
	__m128d x,y;
	y=_mm_sub_pd(p,q);
	x=_mm_div_pd(y,b);
	x=_mm_round_pd(x,0x0);
	x=_mm_mul_pd(x,b);
	x=_mm_sub_pd(x,y);
	return x;	
}
double dist(__m128d p,__m128d q,__m128d b){
	__m128d x,y;
	y=_mm_sub_pd(p,q);
	x=_mm_div_pd(y,b);
	x=_mm_round_pd(x,0x0);
	x=_mm_mul_pd(x,b);
	x=_mm_sub_pd(x,y);
	x=_mm_dp_pd(x,x,0xFF);
	return _mm_cvtsd_f64(x);
}
__m128d _mm_dist_uy(__m128d p,__m128d q,__m128d b,double uy){
	__m128d x,y;
	y=_mm_sub_pd(p,q);
	x=_mm_div_pd(y,b);
	x=_mm_round_pd(x,0x0);
	x=_mm_mul_pd(x,b);
	x=_mm_sub_pd(x,y);
	x[0]+=x[1]*uy;
	return x;	
}
double dist_uy(__m128d p,__m128d q,__m128d b,double uy){
	__m128d x,y;
	y=_mm_sub_pd(p,q);
	x=_mm_div_pd(y,b);
	x=_mm_round_pd(x,0x0);
	x=_mm_mul_pd(x,b);
	x=_mm_sub_pd(x,y);
	x=_mm_dp_pd(x,x,0xFF);
	x[0]+=x[1]*uy;
	return _mm_cvtsd_f64(x);
}
double length2(__m128d x){
	x=_mm_dp_pd(x,x,0xFF);
	return _mm_cvtsd_f64(x);
}
double length(__m128d x){
	x=_mm_dp_pd(x,x,0xFF);
	return sqrt(_mm_cvtsd_f64(x));
}
double dot(__m128d a,__m128d b){
	__m128d c=_mm_dp_pd(a,b,0x31);
	return _mm_cvtsd_f64(c);
}
__m128d normalize(__m128d a){
	double n=length(a);
	return a/n;
}
double *init_dexp(double e){
	int i,j;
	double *dexp=alloc(2*NDEXP*sizeof(double));
	for(i=-NDEXP,j=0;i<NDEXP;i++,j++)*(dexp+j)=exp(i*e);
	return dexp;
}
__m128d rnd11(){
	__m128d b={dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt)};
	return b;
}
__m128d sincosa(double a){
	__m128d b={sin(a),cos(a)};
	return b;
}
__m128d rot22(__m128d a,__m128d b01){
	__m128d b10={-b01[1],b01[0]};
	return _mm_set1_pd(a[1])*b01-_mm_set1_pd(a[0])*b10;
}
__m128d rot2w(__m128d a,double w){
	__m128d b={-sin(w),cos(w)};
	return normalize(rot22(a,b));
}
unsigned int rndn(unsigned int n){
	return (unsigned)(dsfmt_genrand_open_open(&dsfmt)*n);
}
