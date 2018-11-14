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
#include "io.h"
#include "rnd.h"
#include "colors.h"

#if defined SDL
//Graphics
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include "mySDL.h"
#include "draw.h"
#endif

//#include "postscript.h"
extern dsfmt_t dsfmt;
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
void explore(header *t){
	//Explore parameter space
	/////////////////////////
	double rnd;
	//epsilon
	/////////
	rnd=(dsfmt_genrand_open_open(&dsfmt)-0.5)*0.1;
	t->epsilon+=rnd;
	if(t->epsilon<0.0)t->epsilon*=-1;
	//pressure
	/////////
	rnd=(dsfmt_genrand_open_open(&dsfmt)-0.5)*0.1;
	t->pressure+=rnd;
	if(t->pressure<0.0)t->pressure*=-1;
	//chemical potential
	///////////////////
	rnd=(dsfmt_genrand_open_open(&dsfmt)-0.5)*1.0;
	t->specie->mu+=rnd;
}
#if defined SDL

#endif
#if defined SDL
int run(header *t,mySDL *sdl){
#else
int run(header *t){
#endif
	//Main routine -- Running the simulation
	////////////////////////////////////////
	long long int i;
	unsigned int ncycle;
	compound_particle *c;
	open_files(t);
	//Optimization
	//////////////
	int acc_move[2]={0,0};
	int acc_rotate[2]={0,0};
	int acc_volume[2]={0,0};
	int acc_volume_xy[2]={0,0};
	int acc_volume_dxdy[2]={0,0};
	int acc_shape[2]={0,0};
	t->max_displacement=_mm_set1_pd(t->max_displacement[0]);

	int *acc[]={acc_move,acc_rotate,acc_volume,acc_volume_xy,acc_volume_dxdy,acc_shape};
	double *mmax[]={&(t->max_displacement[0]),&t->max_rotation,&t->max_vol,&t->max_xy,&t->max_dxdy,&t->max_uy};
	double frac[6];
	
	int count=0;
	char s_name[1024];
	//Signals
	/////////
	signal(SIGINT,signal_safe_exit_int);
	signal(SIGUSR1,signal_safe_exit);
	time(&t->t1),t->t2=t->t1;

	all_particle_energy_hash(t,&t->energy);
	printf(">>>Running canonical simulation ["RED"%d"RESET"]\n",t->energy);
	//alloc_graph(t);
	print_clusters(t);
	for(i=1;!safe_exit&&i<t->step;i++){
#if defined SDL
		if(t->display){
			handle_event(sdl,&safe_exit);
		}
#endif
		//for(ncycle=0;ncycle<t->nparticle;ncycle++){
		for(ncycle=0;ncycle<10;ncycle++){
			//Monte Carlo cycle
			///////////////////
			c=rnd_compound(t);
			if(0.5<dsfmt_genrand_open_open(&dsfmt)){
				acc_move[mc_move(c,t,&t->energy)]++;	
			}
			if(0.5<dsfmt_genrand_open_open(&dsfmt)){
				acc_rotate[mc_rotate(c,t,&t->energy)]++;
				//acc_rotate[mc_rotate_restricted(c,t,&t->energy)]++;
			}
		}
		//Grand canonical moves
		///////////////////////
		if(0.5>dsfmt_genrand_open_open(&dsfmt)){
				mc_gc(t,&t->energy);
				//mc_gc_restricted(t,&t->energy);
		}
		//NPT
		/////
		if(t->npt){
			acc_volume[mc_npt(t,&t->energy)]++;
		}
		if(t->nptxy){
			acc_volume_xy[mc_npt_xy(t,&t->energy)]++;
		}
		if(t->nptxy){
			acc_volume_dxdy[mc_npt_dxdy(t,&t->energy)]++;
		}
		if(t->shape){
			acc_shape[mc_uy(t,&t->energy)]++; //FIXME
		}
		if(!(i%(t->mod*t->pmod))){
			write_files(t);
			flush_files(t);
			gfrac(acc,frac,6);
			if(t->optimize){
				optimize(mmax,frac,6);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}
#if defined SDL
			if(t->display){
				find_all_clusters(t);
				draw(sdl,t);
				//save_png(s_name,sdl);
			}
#endif
			if(t->verbose){
				time(&t->t2);
				print_log(stdout,t,i,difftime(t->t2,t->t1),frac);
				sprintf(s_name,"%s_%d",t->name,count++);
			}
			if(t->explore){
				explore(t);
			}
		}
	}
#if defined SDL
	if(t->display){
		draw(sdl,t);
		//save_png(t->name,sdl);
		SDL_Quit();
	}
#endif
	write_files(t);
	save_configuration(t->name,t);
	close_files(t);
	//find_all_cycles(t);
	checksum(stdout,t,t->energy);

	return 0;
}
