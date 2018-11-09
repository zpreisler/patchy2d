#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#include "program_gl.h"
#include "mySDL.h"
#include <png.h>
#define MAX_SOURCE_SIZE (0x100000)
void mySDLresize(mySDL *s){
	//Routines for window resizing
	//////////////////////////////
	
	float a;
	float m=s->scale;
	SDL_GetWindowSize(s->window,&s->w,&s->h);
	glViewport(0,0,s->w,s->h);
	a=(float)s->h/(float)s->w;

	//FIXME
	if(a>1.0){
		s->proj_matrix=(float[4]){m,0.0,0.0,m/a};
		s->view_matrix=(float[4]){1,0,0,1};
	}
	else{
		s->proj_matrix=(float[4]){m*a,0.0,0.0,m};
		s->view_matrix=(float[4]){1,0,0,1};
	}

	glUseProgram(s->program[0]);
	glUniformMatrix2fv(s->proj_matrix_loc,1,GL_FALSE,s->proj_matrix);
	glUniformMatrix2fv(s->view_matrix_loc,1,GL_FALSE,s->view_matrix);
	glUniform1f(s->uy_loc,s->uy);

	glUseProgram(s->program[1]);
	glUniformMatrix2fv(s->proj_matrix_loc,1,GL_FALSE,s->proj_matrix);
	glUniformMatrix2fv(s->view_matrix_loc,1,GL_FALSE,s->view_matrix);
	glUniform1f(s->uy_loc,s->uy);

	return;
}
mySDL *mySDLinit(unsigned int w,unsigned int h){
	//Initialize SDL
	////////////////

	//Allocate
	mySDL *s=alloc(sizeof(mySDL));
	s->proj_matrix=alloc(sizeof(float)*4);
	s->view_matrix=alloc(sizeof(float)*4);
	s->box=alloc(sizeof(float)*8);

	//Generate program
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,24);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS,1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,4);

	//Init
	SDL_Init(SDL_INIT_VIDEO);

	//Define window
	s->window=SDL_CreateWindow("SDL",
			SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED,w,h,
			SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);

	//Context
	s->glcontext=SDL_GL_CreateContext(s->window);
	SDL_GetWindowSize(s->window,&s->w,&s->h);

	//glewInit
	s->glew_status=glewInit();

	//Create program -- load shaders
	s->program[0]=create_program("/home/zdenek/Projects/patchy2d/shader.vert",
			"/home/zdenek/Projects/patchy2d/shader.geom","/home/zdenek/Projects/patchy2d/shader.frag");
	s->program[1]=create_program("/home/zdenek/Projects/patchy2d/boundry.vert",
			"/home/zdenek/Projects/patchy2d/boundry.geom","/home/zdenek/Projects/patchy2d/boundry.frag");
	
	//Assign uniforms
	s->proj_matrix_loc=glGetUniformLocation(s->program[0],"proj_matrix");
	s->view_matrix_loc=glGetUniformLocation(s->program[0],"view_matrix");
	s->uy_loc=glGetUniformLocation(s->program[0],"uy");

	//Generate buffers
	glGenVertexArrays(2,s->vao); //Generate vertex arrays
	glGenBuffers(3,s->vbo); //Generate buffer objects

	//Multisample
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_MULTISAMPLE_ARB);
	glDisable(GL_DEPTH_TEST);

	//Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	//Print graphics card info
	glsl_info();
	return s;
}
void mySDLpositions(mySDL *s,float *p,int n){
	//put positions into graphics buffer
	////////////////////////////////////

	//write positions
	glBindVertexArray(s->vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[0]); //Select buffer

	//write to layout 0
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*n*2,p,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);

	return;
}
void mySDLcolors(mySDL *s,float *c,int n){
	//put colors into graphics buffer
	/////////////////////////////////

	//write_colors
	glBindVertexArray(s->vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[1]); //Select buffer

	//write to layout 1
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*4*n,c,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(1);

	return;
}
void mySDLboundary(mySDL *s,float *b){

	//write positions
	glBindVertexArray(s->vao[1]);
	glBindBuffer(GL_ARRAY_BUFFER,s->vbo[3]); //Select buffer

	//write to layout 0
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*8,b,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);

	return;
}
void mySDLdisplay(mySDL *s){
	//Display
	/////////
	
	//Clear
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	//Draw particles
	////////////////
	glUseProgram(s->program[0]);
	glBindVertexArray(s->vao[0]);
	glDrawArrays(GL_POINTS,0,s->n);

	//Draw boundry box
	//////////////////
	glUseProgram(s->program[1]);
	glBindVertexArray(s->vao[1]);
	//glLineWidth(20.0f);
	glDrawArrays(GL_LINE_LOOP,0,4);

	SDL_GL_SwapWindow(s->window);

	return;
}
void mySLDbuffer(mySDL *s){
	//buffer
	////////
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	//Draw
	glUseProgram(s->program[0]);
	glBindVertexArray(s->vao[0]);
	glDrawArrays(GL_POINTS,0,s->n);

	//Draw
	glUseProgram(s->program[1]);
	glBindVertexArray(s->vao[1]);
}
void m128d2float(__m128d *a,float *b,int n){
	int i,j;
	double *c;
	for(i=0;i<n;i++){
		j=2*i;
		c=(double*)(a+i);
		*(b+j)=(float)c[0];
		*(b+j+1)=(float)c[1];
	}
}
void mySDLsetcolor(float *b,float *color,int n){
	int i,j,k;
	for(i=0;i<n;i++){
		j=i*4;
		for(k=0;k<4;k++){
			*(b+j+k)=*(color+k);
		}
	}
}
void save_png(char *name,mySDL *s){
	//Save png
	//////////
	//FIXME
	FILE *f=open_file2(name,".png","w");
	mySLDbuffer(s);
	png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
	png_infop png_info;
	png_info=png_create_info_struct(png_ptr);
	setjmp(png_jmpbuf(png_ptr));
	png_init_io(png_ptr,f);
	png_set_IHDR(png_ptr,png_info,s->w,s->h,8,
			PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);

	unsigned char data[3*s->w*s->h],argb[4*s->w*s->h];
	unsigned char *rows[s->h];
	int i,j,i1,i2;
	glReadPixels(0,0,s->w,s->h,GL_BGRA,GL_UNSIGNED_INT_8_8_8_8,argb);
	for(i=0;i<s->h;++i){
		rows[s->h-i-1]=data+(i*s->h*3);
		for(j=0;j<s->w;++j){
			i1=(i*s->h+j)*3;
			i2=(i*s->h+j)*4;
			data[i1++]=argb[++i2];
			data[i1++]=argb[++i2];
			data[i1++]=argb[++i2];
		}
	}
	png_set_rows(png_ptr,png_info,rows);
	png_write_png(png_ptr,png_info,PNG_TRANSFORM_IDENTITY,NULL);
	png_write_end(png_ptr,png_info);
	png_destroy_write_struct(&png_ptr,NULL);
	close_file(f);
}
