#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "dSFMT.h"
#include "program_gl.h"
#include "mySDL.h"
#include "init.h"
#include "run.h"
dsfmt_t dsfmt;

int main(int argc, char *argv[]){
	//main
	//////
	float color[4]={1.0,0.0,0.0,0.333};
	//Initialize header
	///////////////////
	
	header *t=init_header();
	dsfmt_init_gen_rand(&dsfmt,t->seed);

	//Reading inputs
	////////////////
	
	if(argc==1){
		usage(stdout,t->argz,*argv);
		return 0;
	}
	printf(">>>Initialized default header\n");
	input_files *input=find_configurational_files(argc,argv);
	read_input(argc,argv,input,t);

	//print simulation arguments
	////////////////////////////
	
	dump_args(stdout,t->argz);

	//Graphics
	//////////
	
	mySDL *sdl=NULL;

	//Allocate memory blocks for the graphics
	/////////////////////////////////////////
	
	if(t->display){

		sdl=mySDLinit();
		
		sdl->positions=alloc(sizeof(float)*2*t->nparticle_alloc);
		sdl->colors=alloc(sizeof(float)*4*t->nparticle_alloc);
		sdl->n=t->nparticle;
		sdl->uy=t->uy;
		m128d2float(t->p->q,sdl->positions,sdl->n);
		mySDLsetcolor(sdl->colors,color,sdl->n);

		//Define boundary
		/////////////////
		
		sdl->box=(float[8]){0.0,0.0,t->box[0],0.0,t->box[0],t->box[1],0.0,t->box[1]};
		sdl->scale=1.0/t->box[0];
		
		mySDLresize(sdl);
		mySDLpositions(sdl,sdl->positions,sdl->n);
		mySDLcolors(sdl,sdl->colors,sdl->n);
		mySDLboundary(sdl,sdl->box);
		mySDLdisplay(sdl);

	}

	//Run the simulation
	////////////////////

	run(t,sdl);
	
	return 0;
}
