#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <utils.h>
#include <string.h>
#include <mpi.h>
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
#include "mpi.h"
#include "graph.h"
#include "optimize.h"
#include "canonical.h"
#include "grand_canonical.h"
#include "npt.h"
#include "restricted.h"
#include "cluster.h"
void print_log(FILE *f,header *t,long long int i,double time,double frac[4]){
	double vol=t->box[0]*t->box[1];
	double rho=(double)t->nparticle/vol;
	fprintf(f,UGREEN"LOG\n"RESET
			CYAN"epsilon:"BLUE" %.3f "RESET
			CYAN"pressure:"BLUE" %.3f "RESET
			CYAN"mu:"BLUE" %.3f\n"RESET
			CYAN"step:"BLUE" %Ld "RESET
			CYAN"time:"BLUE" %.0lfs "RESET
			CYAN"mod:"BLUE" %d "RESET
			CYAN"pmod:"BLUE" %d \n"RESET
			YELLOW"energy:"RED" %d "RESET
			YELLOW"U/N:"URED" %.3lf "RESET
			YELLOW"rho:"UGREEN" %.3lf "RESET
			YELLOW"uy:"UBLUE" %.3lf \n"RESET
			GREEN"N compounds:"BLUE" %d\n"RESET
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
			t->epsilon,t->pressure,t->specie->mu,
			i,time,t->mod,t->pmod,
			t->energy,(double)t->energy/t->nparticle,rho,t->uy,
			t->ncompound,
			t->max_displacement[0],t->max_rotation,t->max_vol,t->max_uy,
			frac[0],frac[1],frac[2],frac[3]);
}
void open_files(header *t){
	t->file.fepsilon=open_file2(t->name,".epsilon","w");
	t->file.fmu=open_file2(t->name,".mu","w");
	t->file.fpressure=open_file2(t->name,".pressure","w");
	t->file.fen=open_file2(t->name,".en","w");
	t->file.frho=open_file2(t->name,".rho","w");
	t->file.fvol=open_file2(t->name,".vol","w");
	t->file.fn=open_file2(t->name,".n","w");
	t->file.ftime=open_file2(t->name,".time","w");
	t->file.fstep=open_file2(t->name,".step","w");
	t->file.fncluster=open_file2(t->name,".ncluster","w");
	t->file.fcluster_avg_size=open_file2(t->name,".cluster_avg_size","w");
	t->file.fcluster_max_size=open_file2(t->name,".cluster_max_size","w");
}
void close_files(header *t){
	close_file(t->file.fepsilon);
	close_file(t->file.fmu);
	close_file(t->file.fpressure);
	close_file(t->file.fen);
	close_file(t->file.frho);
	close_file(t->file.fvol);
	close_file(t->file.fn);
	close_file(t->file.ftime);
	close_file(t->file.fstep);
	close_file(t->file.fncluster);
	close_file(t->file.fcluster_avg_size);
	close_file(t->file.fcluster_max_size);
}
void flush_files(header *t){
	fflush(t->file.fepsilon);
	fflush(t->file.fmu);
	fflush(t->file.fpressure);
	fflush(t->file.fen);
	fflush(t->file.frho);
	fflush(t->file.fvol);
	fflush(t->file.fn);
	fflush(t->file.ftime);
	fflush(t->file.fstep);
	fflush(t->file.fncluster);
	fflush(t->file.fcluster_avg_size);
	fflush(t->file.fcluster_max_size);
}
void write_files(header *t){
	double en;
	double rho;
	double vol;
	double n;
	en=(double)t->energy/(double)t->nparticle;
	rho=(double)t->nparticle/(t->box[0]*t->box[1]);
	vol=(t->box[0]*t->box[1]);
	n=t->ncompound;
	//Write
	///////
	uwrite(&t->epsilon,sizeof(double),1,t->file.fepsilon);
	uwrite(&t->specie->mu,sizeof(double),1,t->file.fmu);
	uwrite(&t->pressure,sizeof(double),1,t->file.fpressure);

	uwrite(&en,sizeof(double),1,t->file.fen);
	uwrite(&rho,sizeof(double),1,t->file.frho);
	uwrite(&vol,sizeof(double),1,t->file.fvol);
	//Write cluster data
	////////////////////
	n=t->cluster->ncluster;
	uwrite(&n,sizeof(double),1,t->file.fncluster);
	double x=(double)t->cluster->max_size/(double)t->ncompound;
	//uwrite(&t->cluster->avg_size,sizeof(double),1,t->file.fcluster_avg_size);
	uwrite(&x,sizeof(double),1,t->file.fcluster_avg_size);
	uwrite(&t->cluster->max_size,sizeof(double),1,t->file.fcluster_max_size);
}
