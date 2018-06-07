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
//#include "postscript.h"
/*extern dsfmt_t dsfmt;
static volatile sig_atomic_t safe_exit=0;
void signal_safe_exit(int sig){
	safe_exit=1;
	printf(RESET"["RED"terminating"RESET"] [signal %d]\n",sig);
}
void signal_safe_exit_int(int sig){
	safe_exit=1;
	printf("\n ["RED"terminating"RESET"] [signal %d]\n",sig);
	signal(sig,SIG_DFL);
	return;
}
unsigned int rnd_rank(mpi_world *mpi){
	return (unsigned)(dsfmt_genrand_open_open(&mpi_dsfmt)*mpi->numtasks);
}*/
/*particle *rnd_particle(header *t){
	species *s=t->specie;
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*t->N);
	while(s&&s->N<(rnd+1)){
		rnd-=s->N;
		s=s->next;
	}
	return (particle*)s->p+rnd;
}
particle *rnd_specie(species *s){
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*s->N);
	return (particle*)s->p+rnd;
}*/
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
/*
void mpi_send_receive_rank(mpi_world *mpi,int *send,int *receive){
	*send=rnd_rank(mpi);
	*send=0;
	do{
		*receive=rnd_rank(mpi);
	}while(*send==*receive);
	*receive=3;
}
void assign_block(mpi_world *mpi,mpi_block *b,header *t,int energy,double rnd){
	b->rank=mpi->rank;
	b->energy=energy;
	b->epsilon=t->epsilon;
	b->rnd=rnd;
}*/
/*void print_nvt_log(FILE *f,header *t,long long int i,double time,int energy,double frac[2]){
	double vol=t->box[0]*t->box[1];
	double rho=(double)t->N/vol;
	fprintf(f,UGREEN"NVT\n"RESET
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
			i,time,t->mod,t->pmod,
			energy,(double)energy/t->N,rho,
			t->max_displacement[0],t->max_rotation,
			frac[0],frac[1]);
}*/
/*
int canonical(mpi_world *mpi __attribute__((unused)),header *t){
	long long int i;
	unsigned int ncycle;
	particle *q;
	int energy=0;

	FILE *fen=open_file2(t->name,".en","w");
	double en;
	time_t t1,t2;

	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	t->max_displacement=_mm_set1_pd(t->max_displacement[0]);

	int *acc[]={acc_move,acc_rotate};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation};
	double frac[2];

	signal(SIGINT,signal_safe_exit_int);
	signal(SIGUSR1,signal_safe_exit);
	time(&t1),t2=t1;

	all_particle_energy_hash(t,&energy);
	printf(">>>[%d] Running canonical simulation ["RED"%d"RESET"]\n",mpi->rank,energy);

	alloc_graph(t);

	for(i=1;!safe_exit&&i<t->step;i++){
		for(ncycle=0;ncycle<2*t->N;ncycle++){
			q=rnd_particle(t);
			if(0.5<dsfmt_genrand_open_open(&dsfmt)||!q->npatch){
				acc_move[mc_move(q,t,&energy)]++;	
			}
			else{
				acc_rotate[mc_rotate(q,t,&energy)]++;
			}
		}
		if(!(i%(t->mod*t->pmod))){
			time(&t2);
			en=(double)energy/(double)t->N;
			uwrite(&en,sizeof(double),1,fen);
			gfrac(acc,frac,2);
			if(t->optimize){
				optimize(mmax,frac,2);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
			if(t->verbose){
				print_nvt_log(stdout,t,i,difftime(t2,t1),energy,frac);
			}
		}
	}
	find_all_cycles(t);
	close_file(fen);
	checksum(stdout,t,energy);
	return 0;
}
typedef struct mpi_exchange{
	int label;
	int energy;
	double epsilon;
	double max_displacement;
	double max_rotation;
}mpi_exchange;
void sort_order(int *order,mpi_exchange *a,int n){
	int i,j;
	int k,l;
	for(i=0;i<n;i++){
		order[i]=i;
	}
	for(j=0;j<n;j++){
		for(i=0;i<n-1;i++){
			k=order[i];
			l=order[i+1];
			if((a+k)->epsilon>(a+l)->epsilon){
				SWAP_ORDER(order[i],order[i+1]);
			}
		}
	}
}
void update_exchange(mpi_exchange *a,header *t,int energy,int label){
	a->label=label;
	a->energy=energy;
	a->epsilon=t->epsilon;
	a->max_displacement=t->max_displacement[0];
	a->max_rotation=t->max_rotation;
}
void mpi_create_datatype(MPI_Datatype *mpi_exchange_type){
	MPI_Datatype types[5]={MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	int block_lengths[5]={1,1,1,1,1};
	MPI_Aint offset[5];
	offset[0]=offsetof(mpi_exchange,label);
	offset[1]=offsetof(mpi_exchange,energy);
	offset[2]=offsetof(mpi_exchange,epsilon);
	offset[3]=offsetof(mpi_exchange,max_displacement);
	offset[4]=offsetof(mpi_exchange,max_rotation);
	MPI_Type_create_struct(5,block_lengths,offset,types,mpi_exchange_type);
	MPI_Type_commit(mpi_exchange_type);
}
int mpi_canonical(mpi_world *mpi __attribute__((unused)),header *t){
	long long int i;
	int k;
	unsigned int ncycle;
	particle *q;
	int energy=0;
	double en;
	time_t t1,t2;
	char name[mpi->numtasks][1024];
	char s_name[1024];

	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	t->max_displacement=_mm_set1_pd(t->max_displacement[0]);

	int *acc[]={acc_move,acc_rotate};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation};
	double frac[2];

	FILE *fen=open_file2(t->name,".en","w");

	double rnd;
	double denergy;
	double depsilon;
	int label=mpi->rank;
	int order[mpi->numtasks];

	int count=0;

	mpi_exchange all_exchange[mpi->numtasks];
	MPI_Datatype mpi_exchange_type;
	mpi_create_datatype(&mpi_exchange_type);

	strcpy(name[0],t->name);
	MPI_Allgather(name,128,MPI_CHAR,name,128,MPI_CHAR,MPI_COMM_WORLD);
	update_exchange(all_exchange,t,energy,label);
	MPI_Allgather(all_exchange,1,mpi_exchange_type,all_exchange,1,mpi_exchange_type,MPI_COMM_WORLD);

	sort_order(order,all_exchange,mpi->numtasks);

	signal(SIGINT,signal_safe_exit_int);
	signal(SIGUSR1,signal_safe_exit);
	time(&t1),t2=t1;

	all_particle_energy_hash(t,&energy);
	printf(">>>[%d] Running canonical simulation with replica exchange ["RED"%d"RESET"]\n",mpi->rank,energy);

	for(i=1;!safe_exit&&i<t->step;i++){
		for(ncycle=0;ncycle<2*t->N;ncycle++){
			q=rnd_particle(t);
			if(0.5<dsfmt_genrand_open_open(&dsfmt)||!q->npatch){
				acc_move[mc_move(q,t,&energy)]++;	
			}
			else{
				acc_rotate[mc_rotate(q,t,&energy)]++;
			}
		}
		if(!(i%(t->mod*t->pmod))){
			time(&t2);
			en=(double)energy/(double)t->N;
			uwrite(&en,sizeof(double),1,fen);
			gfrac(acc,frac,2);
			if(t->optimize){
				optimize(mmax,frac,2);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
			if(t->verbose){
				print_nvt_log(stdout,t,i,difftime(t2,t1),energy,frac);
			}
			if(t->snapshot){
				sprintf(s_name,"%s_%d",t->name,count++);
				save_postscript(s_name,t);
				t->steps_passed=i;
				save_configuration(s_name,t);
			}
			update_exchange(all_exchange,t,energy,label);
			MPI_Allgather(all_exchange,1,mpi_exchange_type,all_exchange,1,mpi_exchange_type,MPI_COMM_WORLD);

			rnd=dsfmt_genrand_open_open(&mpi_dsfmt);
			if(rnd<0.5)k=0;
			else k=1;
			for(;k<mpi->numtasks-1;k+=2){
				rnd=dsfmt_genrand_open_open(&mpi_dsfmt);
				depsilon=(all_exchange+order[k])->epsilon-(all_exchange+order[k+1])->epsilon;
				denergy=(all_exchange+order[k])->energy-(all_exchange+order[k+1])->energy;
				if(rnd<exp(-depsilon*denergy)){
					if(mpi->rank==order[k+1]){
						label=(all_exchange+order[k])->label;
						t->epsilon=(all_exchange+order[k])->epsilon;
						t->max_displacement[0]=(all_exchange+order[k])->max_displacement;
						t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
						t->max_rotation=(all_exchange+order[k])->max_rotation;

						fclose(fen);
						strcpy(t->name,name[label]);
						fen=open_file2s(t->name,".en","a");
					}
					else if(mpi->rank==order[k]){
						label=(all_exchange+order[k+1])->label;
						t->epsilon=(all_exchange+order[k+1])->epsilon;
						t->max_displacement[0]=(all_exchange+order[k+1])->max_displacement;
						t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
						t->max_rotation=(all_exchange+order[k+1])->max_rotation;

						fclose(fen);
						strcpy(t->name,name[label]);
						fen=open_file2s(t->name,".en","a");
					}
					SWAP_ORDER(order[k],order[k+1]);
				}
			}
		}
	}
	strcpy(t->name,name[label]);
	close_file(fen);
	checksum(stdout,t,energy);
	t->steps_passed=i;
	return 0;
}
*/
