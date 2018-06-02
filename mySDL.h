#ifndef MYSDL_H
#define MYSDL_H
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
	int w,h;
}mySDL;
extern void mySDLresize(mySDL *s);
extern mySDL *mySDLinit();
extern void mySDLpositions(mySDL *s,double *p,int n);
extern void mySDLcolors(mySDL *s,float *c,int n);
extern void mySDLdisplay(mySDL *s);
#endif
