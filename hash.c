#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <utils.h>
#include "params.h"
#include "hash.h"
void dump_m128i_int32(const int32_t *t){
	printf("%d %d %d %d\n",*t,*(t+1),*(t+2),*(t+3));
}
void dump_m128i_int32_hex(const int32_t *t){
	printf("%x %x %x %x\n",*t,*(t+1),*(t+2),*(t+3));
}
void dump_hash(FILE *f,hash_table *table){
	int i;
	particle *p;
	for(i=0;i<M;i++){
		p=table[i].p;
		while(p){
			fprintf(f,"%d_%p\t",i,p);
			p=p->next;
		}
	}
	fputc('\n',f);
	return;
}
unsigned hash(__m128d x,__m128d h1){
	x=_mm_mul_pd(x,h1);
	return (((unsigned)x[0]<<10)^((unsigned)x[1]));
}
void hash1(header *t){
	__m128d a=_mm_div_pd(t->box,t->hash_cell);
	t->h=_mm_cvttpd_epi32(a);
	a=_mm_cvtepi32_pd(t->h);
	t->h1=_mm_div_pd(a,t->box);
}
void hash_dirs(header *t){
	int i,j;
	for(i=-1,t->ndir=0;i<=1;i++){
		for(j=-1;j<=1;j++){
			*(t->dir+t->ndir++)=_mm_set_epi32(0,0,j,i);
		}
	}
}
void hash_alloc(header *t){
	int i;
	t->table=(hash_table*)alloc(M*sizeof(hash_table));
	for(i=0;i<M;i++)(t->table)[i].list=(particle***)alloc(t->ndir*sizeof(particle**));
}
void hash_lists(header *t){
	int i,x,y;
	int hx,hy,sx,sy;
	unsigned h1,h2;
	__m128i d,shift,mul,cmp;
	__m128i z,h0;
	z=_mm_set1_epi32(0);
	h0=_mm_set_epi32(0,0,-1,-1);
	h0=_mm_add_epi32(t->h,h0);
	//Error handeling
	hx=_mm_extract_epi32(t->h,0);
	hy=_mm_extract_epi32(t->h,1);
	if((hx>1023)|(hy>1023)){
		warn("Box is too large for the hash table");
		t->hash_cell*=2.0;
		hash1(t);
		hash_lists(t);
		return;
	}
	//Assign cell neighbours
	for(x=0;x<=hx;x++){
		for(y=0;y<=hy;y++){
			h1=(x<<10)^(y);
			d=_mm_set_epi32(0,0,y,x);
			for(i=0;i<t->ndir;i++){
				shift=_mm_add_epi32(*(t->dir+i),d);
				cmp=_mm_cmpgt_epi32(shift,h0);
				mul=_mm_and_si128(cmp,t->h);
				shift=_mm_sub_epi32(shift,mul);

				cmp=_mm_cmplt_epi32(shift,z);
				mul=_mm_and_si128(cmp,t->h);
				shift=_mm_add_epi32(shift,mul);

				sx=_mm_extract_epi32(shift,0);
				sy=_mm_extract_epi32(shift,1);

				h2=(sx<<10)^(sy);
				(t->table)[h1].list[i]=&(t->table)[h2].p;
			}
		}
	}
}
void set_hash_cell(header *t){
	species *s=t->specie;
	double d=0;
	while(s){
		d=MAX(d,s->sigma);
		d=MAX(d,s->sigma_well);
		s=s->next;
	}
	t->hash_cell=_mm_set1_pd(d);
}
void set_hash(header *t){
	set_hash_cell(t);
	hash1(t);
	hash_dirs(t);
	hash_alloc(t);
	hash_lists(t);
}
void hash_insert(particle *p,__m128d h1,hash_table *table){
	unsigned h=hash(*(p->q),h1);
	p->next=table[h].p;
	table[h].p=p;
	p->prev=&table[h].p;
	if(p->next)p->next->prev=&p->next;
	p->h=h;
}
void hash_delete(particle *p){
	*p->prev=p->next;
	if(p->next)p->next->prev=p->prev;
}
void hash_reinsert(particle *p,__m128d h1,hash_table *table){
	unsigned h=hash(*(p->q),h1);
	*p->prev=p->next;
	if(p->next)p->next->prev=p->prev;
	p->next=table[h].p;
	table[h].p=p;
	p->prev=&table[h].p;
	if(p->next)p->next->prev=&p->next;
	p->h=h;
}
