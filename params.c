#include <stdio.h>
#include <utils.h>
#include <immintrin.h>
#include <time.h>
#include "zargs.h"
#include "params.h"
/*mpi_world *init_mpi_world(void){
	mpi_world *mpi=(mpi_world*)alloc(sizeof(mpi_world));
	mpi->seed=SEED;
	return mpi;
}*/
header *init_header(void){
	//Initialize header wit hdefualt values
	header *t=(header*)alloc(sizeof(header));
	//Name
	*t->name='a';
	t->specie=NULL;
	t->update=NULL;
	t->nspecies=1;
	//Steps
	t->step=100;
	t->steps_passed=0;
	//HASH
	t->dir=(__m128i*)alloc(DIR_LENGTH*sizeof(__m128i));
	t->box=_mm_set1_pd(8.0);
	t->copy=(__m128d){1.0,1.0};
	t->uy=0.0;
	t->hash_cell=_mm_set1_pd(1.0);
	t->argz.args=NULL;
	//Parameters
	t->epsilon=1.0,t->pressure=1.0;
	t->end_epsilon=0.0;
	t->lambda=1.0,t->lambda_coupling=1.0;
	//
	t->max_displacement=_mm_set1_pd(1.0);
	t->max_rotation=1.0;
	t->max_vol=0.001;
	t->max_xy=0.001;
	t->max_dxdy=0.001;
	t->max_uy=0.1;
	t->max_dsigma=0.01;
	//Wrinting and printing modulus
	t->mod=10,t->pmod=10;
	t->optimize=1;
	t->snapshot=1;
	t->verbose=1;
	t->display=1;
	t->explore=0;
	//Ensembles
	///////////
	t->npt=0;
	t->nptxy=0;
	t->nptdxdy=0;
	t->shape=0;
	//Random seed
	t->seed=1313;
	t->screen_geometry=(__m128d){640,640};
	//assigning arguments
	time(&t->init_time);
	set_args(t);
	return t;
}
