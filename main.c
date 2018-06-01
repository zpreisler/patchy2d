#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#define MAX_SOURCE_SIZE (0x100000)
unsigned int vao[10];
unsigned int vbo[10];

void initialize(){
	unsigned int program;

	program=create_program("shader.vert","shader.geom","shader.frag");

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float v[]={-0.5,-0.5,0.5,-0.5,0.0,0.5};
	float v2[]={-0.5,0.5,0.5,0.5,0.0,-0.5};
	float color[]={
		0.5,1.0,0.5,1.0,
		0.5,1.0,0.75,1.0,
		0.5,1.0,0.35,1.0
	};

	glGenVertexArrays(1,vao);
	glGenBuffers(2,vbo); //Generate buffers

	glBindVertexArray(vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,vbo[0]); //Select buffer
	glBufferData(GL_ARRAY_BUFFER,sizeof(v),v,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER,vbo[1]); //Select buffer
	glBufferData(GL_ARRAY_BUFFER,sizeof(color),color,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(1);
}
void display(){
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glBindVertexArray(vao[0]);
	glDrawArrays(GL_POINTS,0,3);
}

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window *window=SDL_CreateWindow("SDL",
			SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED,640,640,
			SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);

	SDL_GLContext glcontext=SDL_GL_CreateContext(window);
	GLenum glew_status = glewInit();

	glsl_info();
	initialize();

	display();
	SDL_GL_SwapWindow(window);
	SDL_Delay(1000);

	SDL_Quit();
	return 0;
}
