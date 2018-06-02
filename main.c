#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#define MAX_SOURCE_SIZE (0x100000)
typedef struct mySDL{
	SDL_Window *window;
	SDL_GLContext glcontext;
	SDL_Event event;
	GLenum glew_status;
	unsigned int program[2];
	unsigned int vao[1];
	unsigned int vbo[2];
	float *proj_matrix;
	float *view_matrix;
	unsigned int proj_matrix_loc;
	unsigned int view_matrix_loc;
	float box_offset[4];
}mySDL;
mySDL *mySDLInit(){
	//Allocate
	mySDL *s=alloc(sizeof(mySDL));
	s->proj_matrix=alloc(sizeof(float)*4);
	s->view_matrix=alloc(sizeof(float)*4);
	//Generate program
	SDL_Init(SDL_INIT_VIDEO);
	//Define window
	s->window=SDL_CreateWindow("SDL",
			SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED,640,640,
			SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
	s->glcontext=SDL_GL_CreateContext(s->window);
	s->glew_status=glewInit();
	//Create program -- load shaders
	*s->program=create_program("shader.vert","shader.geom","shader.frag");
	//Assign uniforms
	s->proj_matrix_loc=glGetUniformLocation(*s->program,"proj_matrix");
	s->view_matrix_loc=glGetUniformLocation(*s->program,"view_matrix");
	//Generate buffers
	glGenVertexArrays(1,s->vao); //Generate vertex arrays
	glGenBuffers(2,s->vbo); //Generate buffer objects
	//Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//Print graphics card info
	glslInfo();
	return s;
}
void load(mySDL *s){
	//float proj_matrix[4]={1.0,0.0,0.0,1.0};
	//float view_matrix[4]={-1.0,-1.0,2.0,2.0};

	glUniformMatrix2fv(s->proj_matrix_loc,1,GL_FALSE,s->proj_matrix);
	glUniformMatrix2fv(s->view_matrix_loc,1,GL_FALSE,s->view_matrix);

	float v[]={0.0,0.0,1.0,1.0,2.0,2.0};
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
	glDrawArrays(GL_POINTS,0,3);
	SDL_GL_SwapWindow(s->window);
}
int main(int argc, char *argv[]){
	int quit=0;
	mySDL *s=mySDLInit();
	s->proj_matrix=(float[4]){1.0,0.0,0.0,1.0};
	s->view_matrix=(float[4]){0.0,0.0,2.0,2.0};
	int w,h;

	while (!quit){
		SDL_WaitEvent(&s->event);
		switch (s->event.type){
			case SDL_QUIT:
				quit=1;
				break;
			case SDL_WINDOWEVENT:
				if(s->event.window.event==SDL_WINDOWEVENT_RESIZED){
					SDL_GetWindowSize(s->window,&w,&h);
					s->proj_matrix=(float[4]){(float)w/h,0.0,0.0,(float)w/h};
					printf("s w %d h %d\n",w,h);
				}
				break;
		}
		load(s);
		display(s);
	}
	SDL_Quit();

	return 0;
}
