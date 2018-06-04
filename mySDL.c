#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#include "mySDL.h"
#define MAX_SOURCE_SIZE (0x100000)
void mySDLresize(mySDL *s){
	//Routines for window resizing
	float a;
	float m=s->scale;
	SDL_GetWindowSize(s->window,&s->w,&s->h);
	glViewport(0,0,s->w,s->h);
	a=(float)s->h/(float)s->w;
	if(a>1.0){
		s->proj_matrix=(float[4]){m,0.0,0.0,m/a};
		s->view_matrix=(float[4]){s->box[0],s->box[1],s->box[4],s->box[5]};
	}
	else{
		s->proj_matrix=(float[4]){m*a,0.0,0.0,m};
		s->view_matrix=(float[4]){s->box[0],s->box[1],s->box[4],s->box[5]};
	}
	glUseProgram(s->program[0]);
	glUniformMatrix2fv(s->proj_matrix_loc,1,GL_FALSE,s->proj_matrix);
	glUniformMatrix2fv(s->view_matrix_loc,1,GL_FALSE,s->view_matrix);

	glUseProgram(s->program[1]);
	glUniformMatrix2fv(s->proj_matrix_loc,1,GL_FALSE,s->proj_matrix);
	glUniformMatrix2fv(s->view_matrix_loc,1,GL_FALSE,s->view_matrix);
}
mySDL *mySDLinit(){
	//Allocate
	mySDL *s=alloc(sizeof(mySDL));
	s->proj_matrix=alloc(sizeof(float)*4);
	s->view_matrix=alloc(sizeof(float)*4);
	s->box=alloc(sizeof(float)*8);
	//Generate program
	SDL_Init(SDL_INIT_VIDEO);
	//Define window
	s->window=SDL_CreateWindow("SDL",
			SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED,640,480,
			SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
	s->glcontext=SDL_GL_CreateContext(s->window);
	SDL_GetWindowSize(s->window,&s->w,&s->h);
	s->glew_status=glewInit();
	//Create program -- load shaders
	s->program[0]=create_program("shader.vert","shader.geom","shader.frag");
	s->program[1]=create_program("boundry.vert","boundry.geom","boundry.frag");
	//s->program[1]=create_program("boundry.vert",NULL,"boundry.frag");
	//Assign uniforms
	s->proj_matrix_loc=glGetUniformLocation(s->program[0],"proj_matrix");
	s->view_matrix_loc=glGetUniformLocation(s->program[0],"view_matrix");
	//Generate buffers
	glGenVertexArrays(2,s->vao); //Generate vertex arrays
	glGenBuffers(3,s->vbo); //Generate buffer objects
	//Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//Print graphics card info
	glsl_info();
	return s;
}
void mySDLpositions(mySDL *s,float *p,int n){
	//write positions
	glBindVertexArray(s->vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[0]); //Select buffer
	//write to layout 0
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*n,p,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);
}
void mySDLcolors(mySDL *s,float *c,int n){
	//write_colors
	glBindVertexArray(s->vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[1]); //Select buffer
	//write to layout 1
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*n,c,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(1);
}
void mySDLboundary(mySDL *s,float *b,int n){
	//write positions
	glBindVertexArray(s->vao[1]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[3]); //Select buffer
	//write to layout 0
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*n,b,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);
}
void mySDLdisplay(mySDL *s){
	//Draw particles and the boundary
	//Clear
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	//Draw
	glUseProgram(s->program[0]);
	glBindVertexArray(s->vao[0]);
	glDrawArrays(GL_POINTS,0,3);
	//Draw
	glUseProgram(s->program[1]);
	glBindVertexArray(s->vao[1]);
	//glLineWidth(20.0f);
	glDrawArrays(GL_LINE_LOOP,0,4);

	SDL_GL_SwapWindow(s->window);
}
