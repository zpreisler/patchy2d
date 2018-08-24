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

#if defined SDL
//Graphics
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include "mySDL.h"
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
particle *rnd_particle(header *t){
	species *s=t->specie;
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*t->nparticle);
	while(s&&s->nparticle<(rnd+1)){
		rnd-=s->nparticle;
		s=s->next;
	}
	return (particle*)s->p+rnd;
}
compound_particle *rnd_compound(header *t){
	species *s=t->specie;
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*t->ncompound);
	while(s&&s->ncompound<(rnd+1)){
		rnd-=s->ncompound;
		s=s->next;
	}
	return (compound_particle*)s->c+rnd;
}
particle *rnd_specie(species *s){
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*s->nparticle);
	return (particle*)s->p+rnd;
}
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

	//chmical potential
	///////////////////
	
	rnd=(dsfmt_genrand_open_open(&dsfmt)-0.5)*1.0;
	t->specie->mu+=rnd;
}
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

	//Open files
	////////////
	
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
		//SDL sindow and polling events
		///////////////////////////////

		if(t->display){
			SDL_PollEvent(&sdl->event);
			switch(sdl->event.type){
				case SDL_QUIT:
					safe_exit=1;
					break;
				case SDL_WINDOWEVENT:
					if(sdl->event.window.event==SDL_WINDOWEVENT_RESIZED){
						mySDLresize(sdl);
					}
					break;
			}
		}
#endif

		//Monte Carlo cycle
		///////////////////
		
		for(ncycle=0;ncycle<t->nparticle;ncycle++){
			//Translation and Rotation
			c=rnd_compound(t);
			if(0.5<dsfmt_genrand_open_open(&dsfmt)){
				acc_move[mc_move(c,t,&t->energy)]++;	
			}
			if(0.5<dsfmt_genrand_open_open(&dsfmt)){
				//acc_rotate[mc_rotate(c,t,&t->energy)]++;
				acc_rotate[mc_rotate_restricted(c,t,&t->energy)]++;
			}
			if(0.001>dsfmt_genrand_open_open(&dsfmt)){
				mc_gc_restricted(t,&t->energy);
			}
		}

		//Grand canonical moves
		///////////////////////
		
		if(0.5>dsfmt_genrand_open_open(&dsfmt)){
				//mc_gc_restricted(t,&t->energy);
		//		mc_gc(t,&t->energy);
		}
		if(0.5>dsfmt_genrand_open_open(&dsfmt)){
			//acc_volume_xy[mc_npt_xy(t,&t->energy)]++;
			//acc_volume_dxdy[mc_npt_dxdy(t,&energy)]++;
			//acc_shape[mc_uy(t,&energy)]++; //FIXME
			//acc_volume[mc_npt(t,&t->energy)]++;
		}

		//Print
		///////
		
		if(!(i%(t->mod*t->pmod))){
			find_all_clusters(t);

			write_files(t);
			flush_files(t);
			gfrac(acc,frac,6);

			if(t->optimize){
				optimize(mmax,frac,6);
				t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
			}

			if(t->verbose){
				time(&t->t2);
				print_log(stdout,t,i,difftime(t->t2,t->t1),frac);

				char s_name[1024];
				sprintf(s_name,"%s_%d",t->name,count++);
				printf("%s\n",s_name);

#if defined SDL
				if(t->display){
					save_png(s_name,sdl);
				}
#endif
			}

			//explore
			/////////

			if(t->explore){
				explore(t);
			}

#if defined SDL

			//Update screen
			///////////////
			
			if(t->display){
				
				sdl->scale=1.0/t->box[0];
				sdl->n=t->nparticle;
				sdl->uy=t->uy;
				m128d2float(t->p->q,sdl->positions,sdl->n);
				float color[4]={0.0,1.0,0.0,0.333};
				//float color2[4]={0.0,0.0,1.0,0.333};
				mySDLsetcolor(sdl->colors,color,sdl->n);

				for(unsigned int k=0;k<t->specie->ncompound;k++){
					compound_particle *cc=t->specie->c+k;
					double ox,oy;
					ox=(*(cc)->or)[0];
					oy=(*(cc)->or)[1];
					double a=(atan2(ox,oy));
					for(int ii=0;ii<cc->nparticle;ii++){
						*(sdl->colors+ii*4+k*cc->nparticle*4+0)=cos(a+M_PI/3.0)*cos(a+M_PI/3.0)*0.99;
						*(sdl->colors+ii*4+k*cc->nparticle*4+1)=cos(a+2.0*M_PI/3.0)*cos(a+2.0*M_PI/3.0)*0.66;
						*(sdl->colors+ii*4+k*cc->nparticle*4+2)=cos(a)*cos(a)*0.77;
					}
				}

				mySDLpositions(sdl,sdl->positions,sdl->n);
				mySDLcolors(sdl,sdl->colors,sdl->n);

				printf("ncluster [%d] max_cluster [%.0lf] avg_cluster [%.2lf]\n",t->cluster->ncluster,t->cluster->max_size,t->cluster->avg_size);

				set_all_particle_color(t);
				color_all_clusters(t);

				//color_cluster(t->cluster->clusters,color);

				mySDLcolors(sdl,t->particle_colors,t->nparticle);
				//mySDLboundary(s,s->box);
				mySDLresize(sdl);
				mySDLdisplay(sdl);
			}
#endif

		}
	}

	//clean up
	//////////

#if defined SDL
	if(t->display){
		mySDLdisplay(sdl);
		save_png(t->name,sdl);
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
