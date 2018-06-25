#include <stdlib.h>
#include <stdio.h>
#include <utils.h>
#include <math.h>
#include "zargs.h"
#include "alloc.h"
#include "graph.h"
void count_particles(header *t){
	unsigned int cbox;
	species *s=t->specie;	
	t->ncompound=0;
	t->ncompound_alloc=0;
	t->nparticle=0;
	t->nparticle_alloc=0;
	t->npatches=0; //number of patches
	t->npatches_alloc=0; //number of patches allocated
	//copy box multiplier
	cbox=t->copy[0]*t->copy[1];
	//count over all specie (species are ordered in linked list)
	while(s){
		s->ncompound=s->N*cbox; //s->N is number of compounds
		s->nparticle=s->ncompound*s->nppc; //number of particles is number of compounds x number of particles per compound
		s->npatches=s->nparticle*s->npatch; //number of particle x number of patches per particle
		//alloc
		s->ncompound_alloc=max(s->ncompound*M_ALLOC,MIN_ALLOC);
		s->nparticle_alloc=s->ncompound_alloc*s->nppc;
		s->npatches_alloc=s->nparticle_alloc*s->npatch;
		//count
		t->ncompound+=s->ncompound;
		t->nparticle+=s->nparticle;
		t->npatches+=s->npatches;
		//count alloc
		t->ncompound_alloc+=s->ncompound_alloc;
		t->nparticle_alloc+=s->nparticle_alloc;
		t->npatches_alloc+=s->npatches_alloc;
		//move on	
		s=s->next;
	}
}
compound_memory *alloc_cmem(unsigned size){
	compound_memory *cmem=(compound_memory*)alloc(sizeof(compound_memory));	
	cmem->q=(__m128d*)alloc(size*sizeof(__m128d));
	cmem->q_tmp=(__m128d*)alloc(size*sizeof(__m128d));
	cmem->or=(__m128d*)alloc(size*sizeof(__m128d));
	return cmem;
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
void map_cmem(compound_particle *c,compound_memory *cmem,unsigned k){
	c->q=cmem->q+k;
	c->q_tmp=cmem->q_tmp+k;
	c->or=cmem->or+k;
	//id
	c->n=k;
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
	//p->nd_d2=(double*)alloc(NLIST*sizeof(double));
	p->nd_r2=(double*)alloc(NLIST*sizeof(double));
	p->nd_rij=(__m128d*)alloc(NLIST*sizeof(__m128d));
	//id
	p->n=k;
}
void map_smem(patch *s,patch_memory *smem,unsigned k){
	s->angle=smem->angle+k;
	s->q_angle=smem->q_angle+k;
	s->q=smem->q+k;
}
void assign_particle(particle *p,species *s){
	//patch
	p->sigma=s->sigma;
	p->sigma_well=s->sigma_well;
	p->npatch=s->npatch;
	p->patch_width=cos(s->patch_width/180.0*M_PI);
	//Graphs
	p->pass=0;
	p->idx=-1;
	p->npcycles=0;
	p->pcycles=(pcycle**)alloc(sizeof(pcycle*)*s->npatch);
}
particle *alloc_particles(header *t){
	unsigned i,j;
	int k;
	int kcompound=0,kparticle=0,kpatch=0;
	int type=0;
	unsigned int count_patch=0,count_particle=0,count_compound=0;
	compound_particle *c;
	particle *q;
	patch *ss;
	species *s=t->specie;
	count_particles(t);
	//Allocate memory for compounds particles and patches
	t->c=(compound_particle*)alloc(sizeof(compound_particle)*t->ncompound_alloc);
	t->p=(particle*)alloc(sizeof(particle)*t->nparticle_alloc);
	t->s=(patch*)alloc(sizeof(patch)*t->npatches_alloc);
	//Allocate a memory blocks
	compound_memory *cmem=alloc_cmem(t->ncompound_alloc);
	particle_memory *pmem=alloc_pmem(t->nparticle_alloc);
	patch_memory *smem=alloc_smem(t->npatches_alloc);
	//Map the memory
	for(i=0;i<t->ncompound_alloc;i++){
		map_cmem(t->c+i,cmem,i);
	}
	for(i=0;i<t->nparticle_alloc;i++){
		map_pmem(t->p+i,pmem,i);
	}
	for(i=0;i<t->npatches_alloc;i++){
		map_smem(t->s+i,smem,i);
	}
	j=0;
	while(s){
		s->flag=j++; //Flag different compounds
		s->c=t->c+kcompound;
		s->p=t->p+kparticle;
		s->s=t->s+kpatch;
		//
		kcompound+=s->ncompound_alloc;
		kparticle+=s->nparticle_alloc;
		kpatch+=s->npatches_alloc;
		//Loop over compounds
		for(i=0;i<s->ncompound_alloc;i++){
			c=s->c+i; //c -- current compound
			c->specie=s;
			c->p=s->p+count_particle;
			c->s=s->s+count_patch;
			c->id=count_compound;
			c->nparticle=s->nppc;
			count_compound++;
			for(j=0;j<s->nppc;j++){
				q=c->p+j; //q -- current particle
				//linking
				q->patch=s->s+count_patch;
				q->c=c;
				q->specie=s;
				//flags
				q->type=type;
				q->id=count_particle;
				q->flag=s->flag;
				//assing particle
				assign_particle(q,s);
				count_particle++;
				for(k=0;k<s->npatch;k++){
					ss=q->patch+k;
					ss->id=count_patch;
					ss->p=q;
					ss->c=c;
					ss->specie=s;
					count_patch++;
				}
			}
		}
		type++;
		s=s->next;
	}
	return t->p;
}
