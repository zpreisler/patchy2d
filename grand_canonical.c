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
/*extern dsfmt_t mpi_dsfmt;
static volatile sig_atomic_t gc_safe_exit=0;
void signal_safe_exit_gc(int sig){
	gc_safe_exit=1;
	printf(RESET"["RED"terminating"RESET"] [signal %d]\n",sig);
}
void signal_safe_exit_gc_int(int sig){
	gc_safe_exit=1;
	printf("\n ["RED"terminating"RESET"] [signal %d]\n",sig);
	signal(sig,SIG_DFL);
}*/
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
		n=s->N;
		p=(particle*)s->p;
		zeta=exp(s->mu)*vol;
		//Insert particle
		if(rnd<0.5){
			if(s->N<s->Nalloc-1){
				if(!pre_insert_particle(t,p+n,&enn)){
					a=zeta*exp(enn*t->epsilon)/(s->N+1.0);
					//a=zeta*exp(enn*t->epsilon)/(s->N+2.0);
					rnd=dsfmt_genrand_open_open(&dsfmt);
					if(a>rnd){
						post_insert_particle(t,p+n);
						s->N++;
						t->N++;
						t->npatch+=s->npatch;
						*en+=enn;
					}
				}
			}
		}
		//Delete particle
		else{
			if(s->N>1){
				m=dsfmt_genrand_open_open(&dsfmt)*n;
				enn=(p+m)->en_new;
				a=exp(-enn*t->epsilon)*s->N/(zeta);
				//a=exp(-enn*t->epsilon)*(s->N+1)/(zeta);
				rnd=dsfmt_genrand_open_open(&dsfmt);
				if(a>rnd){
					delete_particle(t,p+n-1,p+m);
					s->N--;
					t->N--;
					t->npatch-=s->npatch;
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
	double rho=(double)t->N/vol;
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
			t->N,t->ncycles,t->epsilon,
			i,time,t->mod,t->pmod,
			energy,(double)energy/t->N,rho,
			t->max_displacement[0],t->max_rotation,
			frac[0],frac[1]);
}
int print_n_species(FILE *f,header *t){
	species *s=t->specie;
	while(s){
		fprintf(f,"%d ",s->N);
		s=s->next;
	}
	fprintf(f,"\n");
	return 0;
}
/*
int grand_canonical(mpi_world *mpi __attribute__((unused)),header *t){
	long long int i;
	unsigned int ncycle;
	particle *q;
	int energy;
	int count=0;

	FILE *fen=open_file2(t->name,".en","w");
	FILE *fen_tot=open_file2(t->name,".en_tot","w");
	FILE *fnparticles=open_file2(t->name,".n","w");
	FILE *fncycles=open_file2(t->name,".ncycles","w");
	double en;
	double en_tot;
	double ncycles;
	double nparticles;
	time_t t1,t2;

	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	t->max_displacement=_mm_set1_pd(t->max_displacement[0]);

	int *acc[]={acc_move,acc_rotate};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation};
	double frac[2];

	double depsilon;
	char name[1024];

	signal(SIGINT,signal_safe_exit_gc_int);
	signal(SIGUSR1,signal_safe_exit_gc);
	time(&t1),t2=t1;

	alloc_graph(t);

	depsilon=(t->end_epsilon-t->epsilon)/t->step;

	all_particle_energy_hash(t,&energy);
	printf(">>>[%d] Running canonical simulation ["RED"%d"RESET"]\n",mpi->rank,energy);

	for(i=1;!gc_safe_exit&&i<t->step;i++){
		for(ncycle=0;ncycle<2*t->N;ncycle++){
			q=rnd_particle(t);
			if(0.5<dsfmt_genrand_open_open(&dsfmt)){
				acc_move[mc_move(q,t,&energy)]++;	
			}
			else{
				acc_rotate[mc_rotate(q,t,&energy)]++;
			}
			if(0.001<dsfmt_genrand_open_open(&dsfmt)){
				mc_gc(t,&energy);
			}
		}
		if(!(i%(t->mod*t->pmod))){
			time(&t2);
			en=(double)energy/(double)t->N;
			uwrite(&en,sizeof(double),1,fen);
			en_tot=(double)energy;
			uwrite(&en_tot,sizeof(double),1,fen_tot);
			nparticles=(double)t->N;
			uwrite(&nparticles,sizeof(double),1,fnparticles);
			find_all_cycles(t);
			ncycles=(double)t->ncycles;
			uwrite(&ncycles,sizeof(double),1,fncycles);
			gfrac(acc,frac,2);
			if(t->optimize){
				optimize(mmax,frac,2);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
			if(t->verbose){
				print_muvt_log(stdout,t,i,difftime(t2,t1),energy,frac);
			}
			if(t->snapshot){
				sprintf(name,"%s_%d",t->name,count++);
				save_postscript(name,t);
				save_configuration(name,t);
			}
			fflush(fen);
			fflush(fen_tot);
			fflush(fnparticles);
			fflush(fncycles);
		}
		if(t->end_epsilon!=0.0){
			t->epsilon+=depsilon;
		}
	}
	close_file(fen);
	close_file(fen_tot);
	close_file(fnparticles);
	close_file(fncycles);
	checksum(stdout,t,energy);
	return 0;
}
int grand_canonical_n(mpi_world *mpi __attribute__((unused)),header *t){
	long long int i;
	unsigned int ncycle;
	particle *q;
	int energy;
	int count=0;

	FILE *fen=open_file2(t->name,".en","w");
	FILE *fen_tot=open_file2(t->name,".en_tot","w");
	FILE *fnparticles=open_file2(t->name,".n","w");
	FILE *fncycles=open_file2(t->name,".ncycles","w");
	FILE *fnspecies=open_file2(t->name,".nspecies","w");
	FILE *fcycles_histogram=open_file2(t->name,".cycles_histogram","w");
	double en;
	double en_tot;
	double ncycles;
	double nparticles;
	time_t t1,t2;

	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	t->max_displacement=_mm_set1_pd(t->max_displacement[0]);

	int *acc[]={acc_move,acc_rotate};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation};
	double frac[2];

	double depsilon;
	char name[1024];
	double nn1;

	signal(SIGINT,signal_safe_exit_gc_int);
	signal(SIGUSR1,signal_safe_exit_gc);
	time(&t1),t2=t1;

	//alloc_graph(t);

	depsilon=(t->end_epsilon-t->epsilon)/t->step;

	all_particle_energy_hash(t,&energy);
	printf(">>>[%d] Running Grand-Canonical simulation ["RED"%d"RESET"]\n",mpi->rank,energy);
	//find_all_cycles(t);
	//print_pcycle_length_histogram(fcycles_histogram,t->graph_pcycle,t->ncycles);

	for(i=1;!gc_safe_exit&&i<t->step;i++){
		for(ncycle=0;ncycle<2*t->N;ncycle++){
			q=rnd_particle(t);
			nn1=(double)t->N/(t->N+2);
			if(nn1>dsfmt_genrand_open_open(&dsfmt)){
				if(0.5<dsfmt_genrand_open_open(&dsfmt)){
					acc_move[mc_move(q,t,&energy)]++;	
				}
				else{
					acc_rotate[mc_rotate(q,t,&energy)]++;
				}
			}
			if((1.0-nn1)>dsfmt_genrand_open_open(&dsfmt)){
				mc_gc(t,&energy);
			}
		}
		if(!(i%(t->mod*t->pmod))){
			time(&t2);
			en=(double)energy/(double)t->N;
			uwrite(&en,sizeof(double),1,fen);
			en_tot=(double)energy;
			uwrite(&en_tot,sizeof(double),1,fen_tot);
			nparticles=(double)t->N;
			uwrite(&nparticles,sizeof(double),1,fnparticles);
			//find_all_cycles(t);
			ncycles=(double)t->ncycles;
			uwrite(&ncycles,sizeof(double),1,fncycles);
			print_n_species(fnspecies,t);
			//print_pcycle_length_histogram(fcycles_histogram,t->graph_pcycle,t->ncycles);
			gfrac(acc,frac,2);
			if(t->optimize){
				optimize(mmax,frac,2);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
			if(t->verbose){
				print_muvt_log(stdout,t,i,difftime(t2,t1),energy,frac);
			}
			if(t->snapshot){
				sprintf(name,"%s_%d",t->name,count++);
				save_postscript(name,t);
				t->steps_passed=i;
				save_configuration(name,t);
			}
			fflush(fen);
			fflush(fen_tot);
			fflush(fnparticles);
			fflush(fncycles);
			fflush(fnspecies);
			fflush(fcycles_histogram);
		}
		if(t->end_epsilon!=0.0){
			t->epsilon+=depsilon;
		}
	}
	close_file(fen);
	close_file(fen_tot);
	close_file(fnparticles);
	close_file(fncycles);
	close_file(fnspecies);
	close_file(fcycles_histogram);
	checksum(stdout,t,energy);
	t->steps_passed=i;
	return 0;
}*/
