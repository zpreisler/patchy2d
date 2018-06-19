#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <utils.h>
#include <string.h>
//#include <mpi.h>
#include <time.h>
#include "dSFMT.h"
#include "zargs.h"
#include "params.h"
#include "mm_math.h"
#include "hash.h"
#include "energy.h"
#include "lists.h"
#include "init.h"
#include "patches.h"
//#include "mpi.h"
#include "graph.h"
#include "optimize.h"
#include "canonical.h"
extern dsfmt_t dsfmt;
int mc_move(particle *p,header *t,int *en){
	int eno=0,enn=0;
	int de;
	double rmd;
	__m128d q0=*(p->q);
	__m128d rnd=(__m128d){dsfmt_genrand_open_open(&dsfmt)-0.5,dsfmt_genrand_open_open(&dsfmt)-0.5};
	__m128d dq=t->max_displacement*rnd;
	dq[0]-=dq[1]*t->uy;
	*(p->q)+=dq;
	boundary(p->q,t->box);
	hash_reinsert(p,t->h1,t->table);
	eno=p->en_new;
	list_swap(p);
	if(overlap(p,t)){
		*(p->q)=q0;
		hash_reinsert(p,t->h1,t->table);
		list_swap(p);
		return 1;
	}
	enn=particle_energy_hash2(p);
	de=enn-eno;
	rmd=dsfmt_genrand_open_open(&dsfmt);
	if(rmd<exp(t->epsilon*de)){
		adjust_lists(p);
		*en+=de;
		return 0;
	}
	else{
		*(p->q)=q0;
		hash_reinsert(p,t->h1,t->table);
		list_swap(p);
		return 1;
	}
}
int mc_rotate(particle *p,header *t,int *en){
	int eno=0,enn=0;
	int de;
	double w=t->max_rotation*(dsfmt_genrand_open_open(&dsfmt)-0.5);
	double rnd=dsfmt_genrand_open_open(&dsfmt);
	__m128d or0=*p->or;
	*p->or=rot2w(*p->or,w);
	eno=p->en_new;
	list_swap(p);
	set_patches(p);
	enn=particle_energy_hash(p,t);
	de=enn-eno;
	rnd=dsfmt_genrand_open_open(&dsfmt);
	if(rnd<exp(t->epsilon*de)){
		adjust_lists(p);
		*en+=de;
		return 0;
	}
	else{
		*(p->or)=or0;
		hash_reinsert(p,t->h1,t->table);
		set_patches(p);
		list_swap(p);
		return 1;
	}
}
