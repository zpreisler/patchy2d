#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <utils.h>
#include <string.h>
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
#include "graph.h"
#include "optimize.h"
#include "canonical.h"
extern dsfmt_t dsfmt;
int c_en_old(compound_particle *c){
	int i;
	int en_old=0;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		en_old+=p->en_new;
	}
	return en_old;
}
int c_en_new(compound_particle *c){
	int i;
	int en_new=0;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		//en_new+=particle_energy_hash(p,t);
		//en_new+=particle_energy_hash2m(p,t);
		en_new+=particle_energy_hash2(p);
	}
	return en_new;
}
void c_list_swap(compound_particle *c){
	int i;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		list_swap(p);
	}
}
void c_adjust_lists(compound_particle *c){
	int i;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		adjust_lists(p);
	}
}
int mc_move(compound_particle *c,header *t,int *en){
	int eno=0,enn=0;
	int de;
	double rmd;
	__m128d q0=*(c->q);
	__m128d rnd=(__m128d){dsfmt_genrand_open_open(&dsfmt)-0.5,dsfmt_genrand_open_open(&dsfmt)-0.5};
	__m128d dq=t->max_displacement*rnd;
	dq[0]-=dq[1]*t->uy;
	*(c->q)+=dq;
	reset_particle(c,t); //boundary and hash reinsert
	//boundary(p->q,t->box);
	//hash_reinsert(p,t->h1,t->table);
	//eno=p->en_new;
	//list_swap(p);
	//if(overlap(p,t)){
	eno=c_en_old(c);
	c_list_swap(c);
	if(compound_overlap(c,t)){
		*(c->q)=q0;
		reset_particle(c,t);
		c_list_swap(c);
		//hash_reinsert(p,t->h1,t->table);
		//list_swap(p);
		return 1;
	}
	//enn=particle_energy_hash2(p);
	enn=c_en_new(c);
	de=enn-eno;
	rmd=dsfmt_genrand_open_open(&dsfmt);
	if(rmd<exp(t->epsilon*de)){
		c_adjust_lists(c);
		*en+=de;
		return 0;
	}
	else{
		*(c->q)=q0;
		reset_particle(c,t);
		c_list_swap(c);
		//hash_reinsert(p,t->h1,t->table);
		//list_swap(p);
		return 1;
	}
}
int mc_rotate(compound_particle *c,header *t,int *en){
	int eno=0,enn=0;
	int de;
	double rmd;
	double w=t->max_rotation*(dsfmt_genrand_open_open(&dsfmt)-0.5);
	//double rnd=dsfmt_genrand_open_open(&dsfmt);
	__m128d or0=*c->or;
	*c->or=rot2w(*c->or,w);
	reset_particle(c,t);
	//boundary(p->q,t->box);
	//hash_reinsert(p,t->h1,t->table);
	//eno=p->en_new;
	//list_swap(p);
	//if(overlap(p,t)){
	eno=c_en_old(c);
	c_list_swap(c);
	if(compound_overlap(c,t)){
		*(c->or)=or0;
		reset_particle(c,t);
		c_list_swap(c);
		//hash_reinsert(p,t->h1,t->table);
		//list_swap(p);
		return 1;
	}
	//enn=particle_energy_hash2(p);
	enn=c_en_new(c);
	de=enn-eno;
	rmd=dsfmt_genrand_open_open(&dsfmt);
	if(rmd<exp(t->epsilon*de)){
		c_adjust_lists(c);
		//adjust_lists(p);
		*en+=de;
		return 0;
	}
	else{
		*(c->or)=or0;
		reset_particle(c,t);
		c_list_swap(c);
		//hash_reinsert(p,t->h1,t->table);
		//list_swap(p);
		return 1;
	}
}
/*
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
*/
