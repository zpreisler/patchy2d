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
#include "grand_canonical.h"
extern dsfmt_t dsfmt;
inline double lin_wave(double x,int n){
	double p,w;
	p=M_PI/n;
	w=1.0-(fabs(fmod(x,2*p)/p-1.0));
	return w;
}
double en_bias(__m128d a,__m128d b){
	double cos_theta=dot(a,b);
	if(cos_theta>1.0)cos_theta=1.0;
	return lin_wave(cos_theta,1);
}
int mc_rotate_restricted(compound_particle *c,header *t,int *en){
	int eno=0,enn=0;
	int de;
	int k;
	double rmd;
	double a;
	//double cos_theta;
	double en_bias_new,en_bias_old;
	//double beta=fabs(t->epsilon);
	__m128d or0=*c->or;

	k=dsfmt_genrand_open_open(&dsfmt)*3.0;
	a=(double)k/3.0*M_PI;
	*(c)->or=sincosa(a);

	en_bias_old=en_bias(or0,*c->or_well);
	en_bias_new=en_bias(*c->or,*c->or_well);

	rmd=dsfmt_genrand_open_open(&dsfmt);
	a=exp((en_bias_old-en_bias_new)*t->lambda);
	if(!(rmd<a)){
		*(c->or)=or0;
		return 1;
	}

	reset_particle(c,t);
	eno=c_en_old(c);
	c_list_swap(c);
	if(compound_overlap(c,t)){
		*(c->or)=or0;
		reset_particle(c,t);
		c_list_swap(c);
		return 1;
	}
	enn=c_en_new(c);
	de=enn-eno;
	rmd=dsfmt_genrand_open_open(&dsfmt);
	if(rmd<exp(t->epsilon*de)){
		c_adjust_lists(c);
		*en+=de;
		return 0;
	}
	else{
		*(c->or)=or0;
		reset_particle(c,t);
		c_list_swap(c);
		return 1;
	}
}
int mc_gc_specie_restricted(header *t,species *s,int *en){
	double rnd;
	int enn;
	unsigned int n,m;
	compound_particle *c;
	double vol=t->box[0]*t->box[1];
	double zeta;//=exp(t->mu)*vol;
	double a;
	double en_bias_new;
	//Check it is grand canonical enabled specie
	if(s->grand_canonical){
		rnd=dsfmt_genrand_open_open(&dsfmt);
		n=s->ncompound;
		c=s->c;
		zeta=exp(s->mu)*vol;
		//Insert particle
		if(rnd<0.5){
			if(s->ncompound<s->ncompound_alloc-1){
				if(!pre_insert_compound(t,c+n,&enn)){
					en_bias_new=-en_bias(*(c+n)->or,*(c+n)->or_well)*t->lambda;
					//a=zeta*exp(enn*t->epsilon)/(s->ncompound+1.0);
					a=zeta*exp(enn*t->epsilon+en_bias_new)/(s->ncompound+1.0);
					rnd=dsfmt_genrand_open_open(&dsfmt);
					if(a>rnd){
						*en+=enn;
						post_insert_compound(t,c+n);
						s->N++;
						s->ncompound++;
						s->nparticle+=c->nparticle;
						
						t->ncompound++;
						t->nparticle+=c->nparticle;
					}
				}
			}
		}
		//Delete particle
		else{
			if(s->ncompound>1){
				m=dsfmt_genrand_open_open(&dsfmt)*n;
				c=s->c;
				enn=c_en_old(c+m);
				en_bias_new=-en_bias(*(c+m)->or,*(c+m)->or_well)*t->lambda;
				a=exp(-enn*t->epsilon-en_bias_new)*s->ncompound/(zeta);
				rnd=dsfmt_genrand_open_open(&dsfmt);
				if(a>rnd){
					*en-=enn;
					delete_compound(t,c+n-1,c+m);
					s->N--;
					s->ncompound--;
					s->nparticle-=c->nparticle;

					t->ncompound--;
					t->nparticle-=c->nparticle;
				}
			}
		}
	}
	return 0;
}
int mc_gc_restricted(header *t,int *en){
	species *s=t->specie;
	while(s){
		mc_gc_specie_restricted(t,s,en);
		s=s->next;
	}
	return 0;
}

