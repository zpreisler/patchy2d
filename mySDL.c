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
	float a;
	float m=1.0;
	SDL_GetWindowSize(s->window,&s->w,&s->h);
	glViewport(0,0,s->w,s->h);
	a=(float)s->h/(float)s->w;
	if(a>1.0){
		s->proj_matrix=(float[4]){m,0.0,0.0,m/a};
	}
	else{
		s->proj_matrix=(float[4]){m*a,0.0,0.0,m};
	}
	glUniformMatrix2fv(s->proj_matrix_loc,1,GL_FALSE,s->proj_matrix);
	glUniformMatrix2fv(s->view_matrix_loc,1,GL_FALSE,s->view_matrix);
}
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
			SDL_WINDOWPOS_UNDEFINED,640,480,
			SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
	s->glcontext=SDL_GL_CreateContext(s->window);
	SDL_GetWindowSize(s->window,&s->w,&s->h);
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
