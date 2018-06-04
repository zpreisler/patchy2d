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
#define MAX_SOURCE_SIZE 8096
dsfmt_t dsfmt;
int main(int argc, char *argv[]){
	//Initialize header
	header *t=init_header();
	dsfmt_init_gen_rand(&dsfmt,t->seed);
	//Reading inputs
	if(argc==1){
		usage(stdout,t->argz,*argv);
		return 0;
	}
	printf(">>>Initialized default header\n");
	input_files *input=find_configurational_files(argc,argv);
	read_input(argc,argv,input,t);
	
	dump_args(stdout,t->argz);

	//Graphics
	int quit=0;
	float v[6]={0.5,0.0,1.5,0.0,2.5,0.0};
	float color[12]={
		0.5,1.0,0.5,1.0,
		0.5,1.0,0.75,1.0,
		0.5,1.0,0.35,1.0
	};
	mySDL *s=mySDLinit();
	s->box=(float[8]){0.0,0.0,4.0,0.0,4.0,4.0,0.0,4.0};
	s->scale=s->box[4]/(s->box[4]+1.0);
	mySDLresize(s);
	mySDLpositions(s,v,6);
	mySDLcolors(s,color,12);
	mySDLboundary(s,s->box,8);
	mySDLdisplay(s);
	while(!quit){
		SDL_WaitEvent(&s->event);
		switch(s->event.type){
			case SDL_QUIT:
				quit=1;
				break;
			case SDL_WINDOWEVENT:
				if(s->event.window.event==SDL_WINDOWEVENT_RESIZED){
					mySDLresize(s);
				}
				break;
		}
		mySDLdisplay(s);
	}
	SDL_Quit();
	return 0;
}
