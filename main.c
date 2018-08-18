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
	
	mySDL *s=mySDLinit();

	//Allocate memory blocks for the graphics
	/////////////////////////////////////////
	
	s->positions=alloc(sizeof(float)*2*t->nparticle_alloc);
	s->colors=alloc(sizeof(float)*4*t->nparticle_alloc);
	s->n=t->nparticle;
	s->uy=t->uy;
	m128d2float(t->p->q,s->positions,s->n);
	mySDLsetcolor(s->colors,color,s->n);

	//Define boundary
	/////////////////
	
	s->box=(float[8]){0.0,0.0,t->box[0],0.0,t->box[0],t->box[1],0.0,t->box[1]};
	s->scale=1.0/t->box[0];
	
	mySDLresize(s);
	mySDLpositions(s,s->positions,s->n);
	mySDLcolors(s,s->colors,s->n);
	mySDLboundary(s,s->box);
	mySDLdisplay(s);

	//Run the simulation
	////////////////////

	run(t,s);
	
	/*while(!quit){
		SDL_PollEvent(&s->event);
		switch(s->event.type){
			case SDL_QUIT:
				quit=1;
				break;
			case SDL_WINDOWEVENT:
				if(s->event.window.event==SDL_WINDOWEVENT_RESIZED){
					mySDLresize(s);
					mySDLdisplay(s);
				}
				break;
		}
	}*/
	SDL_Quit();

	return 0;
}
