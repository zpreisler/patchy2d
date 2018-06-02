#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#include "mySDL.h"
#define MAX_SOURCE_SIZE (0x100000)
int main(int argc, char *argv[]){
	int quit=0;
	double v[6]={0.0,0.0,1.0,1.0,0.5,0.5};
	float color[]={
		0.5,1.0,0.5,1.0,
		0.5,1.0,0.75,1.0,
		0.5,1.0,0.35,1.0
	};
	mySDL *s=mySDLinit();
	mySDLpositions(s,v,6);
	mySDLcolors(s,color,12);
	mySDLresize(s);
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
