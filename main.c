#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#include "mySDL.h"
#define MAX_SOURCE_SIZE 8096
int main(int argc, char *argv[]){
	int quit=0;
	float v[6]={0.5,0.0,1.5,0.0,2.5,0.0};
	float color[12]={
		0.5,1.0,0.5,1.0,
		0.5,1.0,0.75,1.0,
		0.5,1.0,0.35,1.0
	};
	float b[8]={0.0,0.0,4.0,0.0,4.0,4.0,0.0,4.0};
	mySDL *s=mySDLinit();
	mySDLresize(s);
	mySDLpositions(s,v,6);
	mySDLcolors(s,color,12);
	mySDLboundary(s,b,8);
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
