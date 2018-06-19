#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <utils.h>
#include <mpi.h>
#include <time.h>
#include "dSFMT.h"
#include "init.h"
#include "zargs.h"
#include "params.h"
#include "mm_math.h"
#include "hash.h"
#include "energy.h"
#include "lists.h"
#include "patches.h"
#include "mpi.h"
#include "graph.h"
#include "optimize.h"
#include "canonical.h"
//#include "postscript.h"
extern dsfmt_t dsfmt;
int pre_insert_particle(header *t,particle *q,int *en){
	double a=2.0*M_PI*dsfmt_genrand_open_open(&dsfmt);
	*(q)->q=t->box*rnd11();	
	*(q)->or=sincosa(a);
	if(!overlap(q,t)){
		set_patches(q);
		*en=particle_energy_hash2(q);
		//particle_energy_hash2(q,t,en);
		return 0;
	}
	return 1;
}
int post_insert_particle(header *t,particle *q){
	hash_insert(q,t->h1,t->table);
	adjust_lists(q);
	return 0;
}
int delete_particle(header *t,particle *p,particle *q){ //p is the last particle and q is the particle selected for deletion
	delete_lists(q);
	hash_delete(q);
	if(p!=q){
		swap_lists(p,q);
		LIST_SWAP(p->new_list,q->new_list);
		hash_delete(p);
		q->en_new=p->en_new;
		*(q)->q=*(p)->q;	
		*(q)->or=*(p)->or;
		set_patches(q);
		hash_insert(q,t->h1,t->table);
	}
	return 0;
}
int mc_gc_specie(header *t,species *s,int *en){
	double rnd=dsfmt_genrand_open_open(&dsfmt);
	int enn;
	unsigned int n,m;
	particle *p;
	double vol=t->box[0]*t->box[1];
	double zeta;//=exp(t->mu)*vol;
	double a;
	if(s->grand_canonical){
		rnd=dsfmt_genrand_open_open(&dsfmt);
		n=s->nparticle;
		p=(particle*)s->p;
		zeta=exp(s->mu)*vol;
		//Insert particle
		if(rnd<0.5){
			if(s->nparticle<s->nparticle_alloc-1){
				if(!pre_insert_particle(t,p+n,&enn)){
					a=zeta*exp(enn*t->epsilon)/(s->nparticle+1.0);
					//a=zeta*exp(enn*t->epsilon)/(s->nparticle+2.0);
					rnd=dsfmt_genrand_open_open(&dsfmt);
					if(a>rnd){
						post_insert_particle(t,p+n);
						s->N++;
						s->nparticle++;
						s->npatches+=s->npatch;
						t->nparticle++;
						t->npatches+=s->npatch;
						*en+=enn;
					}
				}
			}
		}
		//Delete particle
		else{
			if(s->nparticle>1){
				m=dsfmt_genrand_open_open(&dsfmt)*n;
				enn=(p+m)->en_new;
				a=exp(-enn*t->epsilon)*s->nparticle/(zeta);
				//a=exp(-enn*t->epsilon)*(s->nparticle+1)/(zeta);
				rnd=dsfmt_genrand_open_open(&dsfmt);
				if(a>rnd){
					delete_particle(t,p+n-1,p+m);
					s->N--;
					s->nparticle--;
					s->npatches-=s->npatch;
					t->nparticle--;
					t->npatches-=s->npatch;
					*en-=enn;
				}
			}
		}
	}
	return 0;
}
int mc_gc(header *t,int *en){
	species *s=t->specie;
	while(s){
		mc_gc_specie(t,s,en);
		s=s->next;
	}
	return 0;
}
void print_muvt_log(FILE *f,header *t,long long int i,double time,int energy,double frac[2]){
	double vol=t->box[0]*t->box[1];
	double rho=(double)t->nparticle/vol;
	fprintf(f,UGREEN"muVT"RESET" N: %d ncycles: %d epsilon: %lf\n"RESET 
			CYAN"step:"BLUE" %Ld "RESET
			CYAN"time:"BLUE" %.0lfs "RESET
			CYAN"mod:"BLUE" %d "RESET
			CYAN"pmod:"BLUE" %d \n"RESET
			YELLOW"energy:"RED" %d "RESET
			YELLOW"U/N:"URED" %.3lf "RESET
			YELLOW"rho:"GREEN" %.3lf \n"RESET
			UBLUE"parameters\n"RESET
			BLACK"displacement:"GREEN" %.4lf "RESET
			BLACK"rotation:"GREEN" %.4lf \n"RESET
			UBLUE"acceptance\n"RESET
			BLACK"displacement:"PURPLE" %.2lf "RESET
			BLACK"rotation:"PURPLE" %.2lf \n"RESET,
			t->nparticle,t->ncycles,t->epsilon,
			i,time,t->mod,t->pmod,
			energy,(double)energy/t->nparticle,rho,
			t->max_displacement[0],t->max_rotation,
			frac[0],frac[1]);
}
int print_n_species(FILE *f,header *t){
	species *s=t->specie;
	while(s){
		fprintf(f,"%d ",s->nparticle);
		s=s->next;
	}
	fprintf(f,"\n");
	return 0;
}
