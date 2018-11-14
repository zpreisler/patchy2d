#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <utils.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include "dSFMT.h"
#include "zargs.h"
#include "params.h"
#include "mm_math.h"
#include "hash.h"
#include "energy.h"
#include "lists.h"
#include "init.h"
#include "patches.h"
#include "mpi.h"
#include "graph.h"
#include "optimize.h"
#include "canonical.h"
#include "grand_canonical.h"
#include "npt.h"
#include "restricted.h"
#include "cluster.h"
#include "io.h"
#include "rnd.h"
#include "colors.h"

#if defined SDL
//Graphics
#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include "mySDL.h"
#include "draw.h"
void draw(mySDL *sdl,header *t){
	float color[4]={0.0,0.0,0.0,1.000};
	sdl->scale=1.0/t->box[0];
	sdl->n=t->nparticle;
	sdl->uy=t->uy;
	m128d2float(t->p->q,sdl->positions,sdl->n);
	mySDLpositions(sdl,sdl->positions,sdl->n);

	mySDLsetcolor(sdl->colors,color,sdl->n);
	mySDLcolors(sdl,sdl->colors,sdl->n);

	//color_cluster(t->cluster->clusters,color);
	set_all_particle_color(t);
	color_all_clusters(t);
	mySDLcolors(sdl,t->particle_colors,t->nparticle);

	mySDLdisplay(sdl);
}
#endif
