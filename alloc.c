#include <stdlib.h>
#include <stdio.h>
#include <utils.h>
#include <math.h>
#include "zargs.h"
#include "alloc.h"
#include "graph.h"
void count_particles(header *t){
	species *s=t->specie;	
	t->N=0;
	t->Nalloc=0;
	t->npatch=0;
	t->npatch_alloc=0;
	while(s){
		s->Nalloc=MAX(s->N*M_ALLOC*t->copy[0]*t->copy[1],MIN_ALLOC);
		t->N+=s->N*t->copy[0]*t->copy[1];
		t->Nalloc+=s->Nalloc;
		t->npatch+=s->N*s->npatch;
		t->npatch_alloc+=s->Nalloc*s->npatch;
		s=s->next;
	}
}
particle_memory *alloc_pmem(unsigned size){
	particle_memory *pmem=(particle_memory*)alloc(sizeof(particle_memory));
	pmem->q=(__m128d*)alloc(size*sizeof(__m128d));
	pmem->q_tmp=(__m128d*)alloc(size*sizeof(__m128d));
	pmem->q_track=(__m128d*)alloc(size*sizeof(__m128d));
	pmem->q_well=(__m128d*)alloc(size*sizeof(__m128d));
	pmem->qp_rij=(__m128d*)alloc(size*sizeof(__m128d));
	pmem->or=(__m128d*)alloc(size*sizeof(__m128d));
	pmem->or_well=(__m128d*)alloc(size*sizeof(__m128d));
	return pmem;
}
patch_memory *alloc_smem(unsigned size){
	patch_memory *smem=(patch_memory*)alloc(sizeof(patch_memory));
	smem->angle=(double*)alloc(size*sizeof(double));
	smem->q_angle=(__m128d*)alloc(size*sizeof(__m128d));
	smem->q=(__m128d*)alloc(size*sizeof(__m128d));
	return smem;
}
void map_pmem(particle *p,particle_memory *pmem,unsigned k){
	p->q=pmem->q+k;
	p->q_tmp=pmem->q_tmp+k;
	p->q_track=pmem->q_track+k;
	p->q_well=pmem->q_well+k;
	p->qp_rij=pmem->qp_rij+k;
	p->or=pmem->or+k;
	p->or_well=pmem->or_well+k;
	//allocate lists
	p->new_list=(particle**)alloc(NLIST*sizeof(particle*));
	p->old_list=(particle**)alloc(NLIST*sizeof(particle*));
	p->nd_list=(particle**)alloc(NLIST*sizeof(particle*));
	p->nd_d2=(double*)alloc(NLIST*sizeof(double));
	//id
	p->n=k;
}
void map_smem(patch *s,patch_memory *smem,unsigned k){
	s->angle=smem->angle+k;
	s->q_angle=smem->q_angle+k;
	s->q=smem->q+k;
}
void assign_particle(particle *p,species *s){
	p->sigma=s->sigma;
	p->sigma_well=s->sigma_well;
	p->npatch=s->npatch;
	p->patch_width=cos(s->patch_width/180.0*M_PI);
	p->pass=0;
	p->idx=-1;
	p->npcycles=0;
	p->flag=s->flag;
	p->pcycles=(pcycle**)alloc(sizeof(pcycle*)*s->npatch);
	p->specie=s;
}
void assign_patch(patch *patch,int id){
	patch->id=id;
	return;
}
particle *alloc_particles(header *t){
	unsigned i,k=0,l=0;
	int j;
	int type=0;
	particle *q;
	species *s=t->specie;
	count_particles(t);
	t->p=(particle*)alloc(sizeof(particle)*t->Nalloc);
	t->s=(patch*)alloc(sizeof(patch)*t->npatch_alloc);
	particle_memory *pmem=alloc_pmem(t->Nalloc);
	patch_memory *smem=alloc_smem(t->npatch_alloc);
	for(i=0;i<t->Nalloc;i++){
		map_pmem(t->p+i,pmem,i);
	}
	for(i=0;i<t->npatch_alloc;i++){
		map_smem(t->s+i,smem,i);
	}
	j=0;
	while(s){
		s->flag=j++;
		s->p=t->p+k;
		k+=s->Nalloc;
		for(i=0;i<s->Nalloc;i++){
			q=s->p+i;
			assign_particle(q,s);
			q->patch=t->s+l;
			l+=s->npatch;
			q->type=type;
			for(j=0;j<s->npatch;j++){
				assign_patch(q->patch+j,l+j);
			}
		}
		type++;
		s=s->next;
	}
	return t->p;
}
