#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include <utils.h>
#include <string.h>
#include <mpi.h>
#include <signal.h>
#include <time.h>
#include "dSFMT.h"
#include "params.h"
#include "zargs.h"
#include "hash.h"
#include "mm_math.h"
#include "patches.h"
#include "energy.h"
#include "lists.h"
#include "alloc.h"
#include "graph.h"
#include "canonical.h"
#include "optimize.h"
#include "init.h"
#include "postscript.h"
extern dsfmt_t dsfmt;
extern dsfmt_t mpi_dsfmt;
static volatile sig_atomic_t safe_exit=0;
void signal_safe_exit_npt(int sig){
	safe_exit=1;
	printf(RESET"["RED"terminating NPT"RESET"] [signal %d]\n",sig);
}
void signal_safe_exit_int_npt(int sig){
	safe_exit=1;
	printf("\n ["RED"terminating"RESET"] [signal %d]\n",sig);
	signal(sig,SIG_DFL);
	return;
}
int mc_npt(header *t,int *en){
	unsigned int i;
	int enn,eno,de;
	double vol,vol_new,xnew;
	double rnd=dsfmt_genrand_open_open(&dsfmt)-0.5;
	double dv=t->max_vol*rnd*t->N;
	double acc;
	double bp;
	species *s;
	particle *q;
	vol=t->box[0]*t->box[1];
	vol_new=vol+dv;
	xnew=sqrt(vol_new/vol);
	__m128d box=t->box;
	eno=*en;
	t->box*=xnew;
	hash1(t);
	hash_lists(t); //FIXME but only if number of the cells changes
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			*(q)->q_tmp=*(q)->q;
			*(q)->q*=xnew;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	int ennb=0;
	ennb=all_particle_energy_hash(t,&enn);
	de=enn-eno;
	bp=t->pressure*fabs(t->epsilon); //betap
	if(ennb!=-1){
		acc=de*t->epsilon-bp*dv+t->N*log(vol_new/vol);
		rnd=dsfmt_genrand_open_open(&dsfmt);
		if(rnd<exp(acc)){
			*en=enn;
			return 0;
		}
	}
	t->box=box;
	hash1(t);
	hash_lists(t);
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			*(q)->q=*(q)->q_tmp;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	return 1;
}
int mc_npt_xy(header *t,int *en){
	unsigned int i;
	int enn,eno,de;
	double vol,vol_new;//,xnew;
	double rnd=dsfmt_genrand_open_open(&dsfmt)-0.5;
	double dx=t->max_vol*rnd;
	double dv;
	double acc;
	double bp;
	unsigned int d=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*2);
	species *s;
	particle *q;
	__m128d box=t->box;
	__m128d f;
	eno=*en;
	vol=t->box[0]*t->box[1];
	if(t->box[d]+dx<3*t->hash_cell[d]){
		return 1;
	}
	t->box[d]+=dx;
	f=t->box/box;
	vol_new=t->box[0]*t->box[1];
	dv=vol_new-vol;
	hash1(t);
	hash_lists(t); //FIXME but only if number of the cells changes
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			*(q)->q_tmp=*(q)->q;
			*(q)->q*=f;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	int ennb=0;
	ennb=all_particle_energy_hash(t,&enn);
	de=enn-eno;
	bp=t->pressure*fabs(t->epsilon); //betap
	if(ennb!=-1){
		acc=de*t->epsilon-bp*dv+t->N*log(vol_new/vol);
		rnd=dsfmt_genrand_open_open(&dsfmt);
		if(rnd<exp(acc)){
			*en=enn;
			return 0;
		}
	}
	t->box=box;
	hash1(t);
	hash_lists(t);
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			*(q)->q=*(q)->q_tmp;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	return 1;
}
int mc_npt_dxdy(header *t,int *en){
	unsigned int i;
	int enn,eno,de;
	double vol,vol_new;//,xnew;
	double rnd=dsfmt_genrand_open_open(&dsfmt)-0.5;
	double dx=t->max_vol*rnd;
	double dv;
	double acc;
	double bp;
	//unsigned int d=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*2);
	species *s;
	particle *q;
	__m128d box=t->box;
	__m128d f;
	eno=*en;
	vol=t->box[0]*t->box[1];
	t->box[0]+=dx;
	t->box[1]=vol/t->box[0];
	f=t->box/box;
	vol_new=t->box[0]*t->box[1];
	dv=vol_new-vol;
	hash1(t);
	hash_lists(t); //FIXME but only if number of the cells changes
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			*(q)->q_tmp=*(q)->q;
			*(q)->q*=f;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	int ennb=0;
	ennb=all_particle_energy_hash(t,&enn);
	de=enn-eno;
	bp=t->pressure*fabs(t->epsilon); //betap
	if(ennb!=-1){
		acc=de*t->epsilon-bp*dv+t->N*log(vol_new/vol);
		rnd=dsfmt_genrand_open_open(&dsfmt);
		if(rnd<exp(acc)){
			*en=enn;
			return 0;
		}
	}
	t->box=box;
	hash1(t);
	hash_lists(t);
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			*(q)->q=*(q)->q_tmp;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	return 1;
}

int mc_uy(header *t,int *en){
	unsigned int i;
	int enn,eno,de;
	double rnd=dsfmt_genrand_open_open(&dsfmt)-0.5;
	double dy=t->max_uy*rnd;
	double acc;
	double uy_old=t->uy;
	species *s;
	particle *q;
	eno=*en;
	t->uy+=dy;
	if(fabs(t->uy)>0.7){
		t->uy=uy_old;
		return 1;
	}
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	int ennb=0;
	ennb=all_particle_energy_hash(t,&enn);
	de=enn-eno;
	if(ennb!=-1){
		acc=de*t->epsilon;
		rnd=dsfmt_genrand_open_open(&dsfmt);
		if(rnd<exp(acc)){
			*en=enn;
			return 0;
		}
	}
	t->uy=uy_old;
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			q=(particle*)s->p+i;
			list_swap(q);
			hash_reinsert(q,t->h1,t->table);
		}
		s=s->next;
	}
	return 1;
}
void print_npt_log(FILE *f,header *t,long long int i,double time,int energy,double frac[4]){
	double vol=t->box[0]*t->box[1];
	double rho=(double)t->N/vol;
	fprintf(f,UGREEN"NPT\n"RESET
			CYAN"step:"BLUE" %Ld "RESET
			CYAN"time:"BLUE" %.0lfs "RESET
			CYAN"mod:"BLUE" %d "RESET
			CYAN"pmod:"BLUE" %d \n"RESET
			YELLOW"energy:"RED" %d "RESET
			YELLOW"U/N:"URED" %.3lf "RESET
			YELLOW"rho:"UGREEN" %.3lf "RESET
			YELLOW"uy:"UBLUE" %.3lf \n"RESET
			UBLUE"parameters\n"RESET
			BLACK"displacement:"GREEN" %.4lf "RESET
			BLACK"rotation:"GREEN" %.4lf "RESET
			BLACK"volume:"GREEN" %.4lf "RESET
			BLACK"shape:"GREEN" %.4lf \n"RESET
			UBLUE"acceptance\n"RESET
			BLACK"displacement:"PURPLE" %.2lf "RESET
			BLACK"rotation:"PURPLE" %.2lf "RESET
			BLACK"volume:"PURPLE" %.2lf "RESET
			BLACK"shape:"PURPLE" %.2lf \n"RESET,
			i,time,t->mod,t->pmod,
			energy,(double)energy/t->N,rho,t->uy,
			t->max_displacement[0],t->max_rotation,t->max_vol,t->max_uy,
			frac[0],frac[1],frac[2],frac[3]);
}
int npt(mpi_world *mpi __attribute__((unused)),header *t){
	long long int i;
	unsigned int ncycle;
	particle *q;
	int energy;

	FILE *fen=open_file2(t->name,".en","w");
	FILE *frho=open_file2(t->name,".rho","w");
	FILE *fvol=open_file2(t->name,".vol","w");
	double en;
	double rho;
	double vol;
	time_t t1,t2;

	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	int acc_volume[2]={0,0};
	int acc_shape[2]={0,0};
	int *acc[]={acc_move,acc_rotate,acc_volume,acc_shape};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation,&t->max_vol,&t->max_uy};
	double frac[4];

	signal(SIGINT,signal_safe_exit_int_npt);
	signal(SIGUSR1,signal_safe_exit_npt);
	time(&t1),t2=t1;

	all_particle_energy_hash(t,&energy);
	printf(">>>[%d] Running NPT simulation\n",mpi->rank);

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
#ifdef NPT_SOLID
		acc_volume[mc_npt_xy(t,&energy)]++;
		acc_shape[mc_uy(t,&energy)]++;
#else
		acc_volume[mc_npt(t,&energy)]++;
#endif
		if(!(i%(t->mod*t->pmod))){
			time(&t2);
			en=(double)energy/(double)t->N;
			uwrite(&en,sizeof(double),1,fen);
			rho=(double)t->N/(t->box[0]*t->box[1]);
			uwrite(&rho,sizeof(double),1,frho);
			vol=(t->box[0]*t->box[1]);
			uwrite(&vol,sizeof(double),1,fvol);
			gfrac(acc,frac,4);
			if(t->optimize){
				optimize(mmax,frac,4);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
			if(t->verbose){
				print_npt_log(stdout,t,i,difftime(t2,t1),energy,frac);
			}
		}
	}
	close_file(fen);
	close_file(frho);
	close_file(fvol);
	checksum(stdout,t,energy);
	return 0;
}
typedef struct mpi_exchange_npt{
	int label;
	int energy;
	double epsilon;
	double pressure;
	double volume;
	double max_displacement;
	double max_rotation;
	double max_vol;
	double max_uy;
}mpi_exchange_npt;
void sort_order_npt(int *order,mpi_exchange_npt *a,int n){
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
void update_exchange_npt(mpi_exchange_npt *a,header *t,int energy,int label){
	a->label=label;
	a->energy=energy;
	a->epsilon=t->epsilon;
	a->pressure=t->pressure;
	a->volume=t->box[0]*t->box[1];
	a->max_displacement=t->max_displacement[0];
	a->max_rotation=t->max_rotation;
	a->max_vol=t->max_vol;
	a->max_uy=t->max_uy;
}
void mpi_create_datatype_npt(MPI_Datatype *mpi_exchange_type){
	MPI_Datatype types[9]={MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	int block_lengths[9]={1,1,1,1,1,1,1,1,1};
	MPI_Aint offset[9];
	offset[0]=offsetof(mpi_exchange_npt,label);
	offset[1]=offsetof(mpi_exchange_npt,energy);
	offset[2]=offsetof(mpi_exchange_npt,epsilon);
	offset[3]=offsetof(mpi_exchange_npt,pressure);
	offset[4]=offsetof(mpi_exchange_npt,volume);
	offset[5]=offsetof(mpi_exchange_npt,max_displacement);
	offset[6]=offsetof(mpi_exchange_npt,max_rotation);
	offset[7]=offsetof(mpi_exchange_npt,max_vol);
	offset[8]=offsetof(mpi_exchange_npt,max_uy);
	MPI_Type_create_struct(9,block_lengths,offset,types,mpi_exchange_type);
	MPI_Type_commit(mpi_exchange_type);
}
int mpi_npt(mpi_world *mpi __attribute__((unused)),header *t){
	long long int i;
	unsigned int ncycle;
	particle *q;
	int k;
	int energy;
	char name[mpi->numtasks][128];
	int count=0;
	char snapshot_name[1024];

	FILE *fen=open_file2(t->name,".en","w");
	FILE *frho=open_file2(t->name,".rho","w");
	FILE *fvol=open_file2(t->name,".vol","w");
	double en;
	double rho;
	double vol;
	time_t t1,t2;

	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	int acc_volume[2]={0,0};
	int acc_shape[2]={0,0};
	int *acc[]={acc_move,acc_rotate,acc_volume,acc_shape};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation,&t->max_vol,&t->max_uy};
	double frac[4];

	double rnd;
	double denergy;
	double depsilon;
	double dbeta_pressure;
	double dvol;
	int label=mpi->rank;
	int order[mpi->numtasks];

	mpi_exchange_npt all_exchange[mpi->numtasks];
	MPI_Datatype mpi_exchange_type;
	mpi_create_datatype_npt(&mpi_exchange_type);

	strcpy(name[0],t->name);
	MPI_Allgather(name,128,MPI_CHAR,name,128,MPI_CHAR,MPI_COMM_WORLD);
	update_exchange_npt(all_exchange,t,energy,label);
	MPI_Allgather(all_exchange,1,mpi_exchange_type,all_exchange,1,mpi_exchange_type,MPI_COMM_WORLD);

	sort_order_npt(order,all_exchange,mpi->numtasks);

	signal(SIGINT,signal_safe_exit_int_npt);
	signal(SIGUSR1,signal_safe_exit_npt);
	time(&t1),t2=t1;

	all_particle_energy_hash(t,&energy);
	printf(">>>[%d] Running NPT simulation\n",mpi->rank);

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
#ifdef NPT_SOLID
		acc_volume[mc_npt_xy(t,&energy)]++;
		acc_shape[mc_uy(t,&energy)]++;
#else
		acc_volume[mc_npt(t,&energy)]++;
#endif
		if(!(i%(t->mod*t->pmod))){
			time(&t2);
			en=(double)energy/(double)t->N;
			uwrite(&en,sizeof(double),1,fen);
			rho=(double)t->N/(t->box[0]*t->box[1]);
			uwrite(&rho,sizeof(double),1,frho);
			vol=(t->box[0]*t->box[1]);
			uwrite(&vol,sizeof(double),1,fvol);
			gfrac(acc,frac,4);
			if(t->optimize){
				optimize(mmax,frac,4);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
			if(t->verbose){
				print_npt_log(stdout,t,i,difftime(t2,t1),energy,frac);
			}
			if(t->snapshot){
				sprintf(snapshot_name,"%s_%d",t->name,count++);
				save_postscript(snapshot_name,t);
				t->steps_passed=i;
				save_configuration(snapshot_name,t);
			}
			//MPI part
			update_exchange_npt(all_exchange,t,energy,label);
			MPI_Allgather(all_exchange,1,mpi_exchange_type,all_exchange,1,mpi_exchange_type,MPI_COMM_WORLD);

			rnd=dsfmt_genrand_open_open(&mpi_dsfmt);
			if(rnd<0.5)k=0;
			else k=1;
			for(;k<mpi->numtasks-1;k+=2){
				rnd=dsfmt_genrand_open_open(&mpi_dsfmt);
				depsilon=(all_exchange+order[k])->epsilon-(all_exchange+order[k+1])->epsilon;
				denergy=(all_exchange+order[k])->energy-(all_exchange+order[k+1])->energy;
				dbeta_pressure=(all_exchange+order[k])->pressure*(all_exchange+order[k])->epsilon
					-(all_exchange+order[k+1])->pressure*(all_exchange+order[k+1])->epsilon;
				dvol=(all_exchange+order[k])->volume-(all_exchange+order[k+1])->volume;
				if(rnd<exp(-depsilon*denergy+dbeta_pressure*dvol)){
					if(mpi->rank==order[k+1]){
						label=(all_exchange+order[k])->label;
						t->epsilon=(all_exchange+order[k])->epsilon;
						t->pressure=(all_exchange+order[k])->pressure;
						t->max_displacement[0]=(all_exchange+order[k])->max_displacement;
						t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
						t->max_rotation=(all_exchange+order[k])->max_rotation;
						t->max_vol=(all_exchange+order[k])->max_vol;
						t->max_uy=(all_exchange+order[k])->max_uy;

						close_file(fen);
						close_file(frho);
						close_file(fvol);
						strcpy(t->name,name[label]);
						fen=open_file2s(t->name,".en","a");
						frho=open_file2s(t->name,".rho","a");
						fvol=open_file2s(t->name,".vol","a");
					}
					else if(mpi->rank==order[k]){
						label=(all_exchange+order[k+1])->label;
						t->epsilon=(all_exchange+order[k+1])->epsilon;
						t->pressure=(all_exchange+order[k+1])->pressure;
						t->max_displacement[0]=(all_exchange+order[k+1])->max_displacement;
						t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
						t->max_rotation=(all_exchange+order[k+1])->max_rotation;
						t->max_vol=(all_exchange+order[k+1])->max_vol;
						t->max_uy=(all_exchange+order[k+1])->max_uy;

						close_file(fen);
						close_file(frho);
						close_file(fvol);
						strcpy(t->name,name[label]);
						fen=open_file2s(t->name,".en","a");
						frho=open_file2s(t->name,".rho","a");
						fvol=open_file2s(t->name,".vol","a");
					}
					SWAP_ORDER(order[k],order[k+1]);
				}
			}
			//
		}
	}
	close_file(fen);
	close_file(frho);
	close_file(fvol);
	checksum(stdout,t,energy);
	t->steps_passed=i;
	return 0;
}

