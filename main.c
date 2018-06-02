#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#include "mySDL.h"
#define MAX_SOURCE_SIZE (0x100000)
void load(mySDL *s){
	float v[]={0.0,0.0,1.0,1.0,0.5,0.5};
	float color[]={
		0.5,1.0,0.5,1.0,
		0.5,1.0,0.75,1.0,
		0.5,1.0,0.35,1.0
	};
	//write positions
	glBindVertexArray(s->vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[0]); //Select buffer
	glBufferData(GL_ARRAY_BUFFER,sizeof(v),v,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);
	//write_colors
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[1]); //Select buffer
	glBufferData(GL_ARRAY_BUFFER,sizeof(color),color,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(1);
}
void display(mySDL *s){
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glDrawArrays(GL_POINTS,0,1);
	SDL_GL_SwapWindow(s->window);
}
int main(int argc, char *argv[]){
	int quit=0;
	mySDL *s=mySDLInit();
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
		load(s);
		display(s);
	}
	SDL_Quit();

	return 0;
}
