#ifndef MYSDL_H
#define MYSDL_H
#include <signal.h>
typedef struct mySDL{
	SDL_Window *window;
	SDL_GLContext glcontext;
	SDL_Event event;
	GLenum glew_status;
	unsigned int program[2];
	unsigned int vao[2];
	unsigned int vbo[3];
	float *proj_matrix;
	float *view_matrix;
	unsigned int proj_matrix_loc;
	unsigned int view_matrix_loc;
	unsigned int uy_loc;
	float uy;
	float *box;
	float scale;
	int w,h;
	//Objects
	float *positions;
	float *colors;
	int n;
}mySDL;
extern void mySDLresize(mySDL *s);
extern mySDL *mySDLinit();
extern void mySDLpositions(mySDL *s,float *p,int n);
extern void mySDLcolors(mySDL *s,float *c,int n);
extern void mySDLboundary(mySDL *s,float *b);
extern void mySDLdisplay(mySDL *s);
extern void m128d2float(__m128d *a,float *b,int n);
extern void mySDLsetcolor(float *b,float *color,int n);

extern void save_png(char *name,mySDL *s);

extern void handle_event(mySDL *sdl,volatile sig_atomic_t *safe_exit);
#endif
