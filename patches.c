#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include "mm_math.h"
#include "zargs.h"
void init_patches(header *t){
	unsigned int i;
	int j;
	patch *p;
	particle *q;
	species *s=t->specie;
	double w;
	__m128d a={0.0,1.0};
	while(s){
		if(!strcmp(s->patch_type,"symmetric")){
			w=2.0*M_PI/s->npatch;
			for(i=0;i<s->Nalloc;i++){
				q=s->p+i;
				for(j=0;j<q->npatch;j++){
					p=q->patch+j;				
					*(p)->q_angle=rot2w(a,w*j);
					*(p)->angle=w*j;
				}
			}
		}
		else if(!strcmp(s->patch_type,"asymmetric")){
			//w=(2.0*M_PI*(360.0-s->patch_angles[0])/360.0)/(s->npatch-1);
			w=(2.0*M_PI*(360.0-s->patch_angle)/360.0)/(s->npatch-1);
			for(i=0;i<s->Nalloc;i++){
				q=s->p+i;
				for(j=0;j<q->npatch;j++){
					p=q->patch+j;				
					*(p)->q_angle=rot2w(a,w*j);
					*(p)->angle=w*j;
				}
			}
		}
		else if(!strcmp(s->patch_type,"one_three")){
			double wj_[s->npatch];
			double wj[s->npatch];
			if(s->npatch>3){
				w=(2.0*M_PI*(360.0-s->angles[0]-s->angles[1])/360.0)/(s->npatch-2);
				for(j=0;j<s->npatch;j++){
					wj_[j]=w;
					wj[j]=0.0;
				}
				wj_[1]=2.0*M_PI*s->angles[0]/360.0;
				wj_[3]=2.0*M_PI*s->angles[1]/360.0;
				for(j=1;j<s->npatch;j++){
					wj[j]=wj[j-1]+wj_[j];
					//printf("angle[%d] %lf %lf\n",j,wj[j],wj[j]/(2.0*M_PI)*360);
				}
			}
			else w=2.0*M_PI/s->npatch;
			for(i=0;i<s->Nalloc;i++){
				q=s->p+i;
				for(j=0;j<q->npatch;j++){
					p=q->patch+j;				
					*(p)->q_angle=rot2w(a,wj[j]);
					*(p)->angle=wj[j];
				}
			}
		}
		s=s->next;
	}
}
void set_patches(particle *p){
	int i;
	patch *s=p->patch;
	for(i=0;i<p->npatch;i++){
		*(s+i)->q=rot22(*(s+i)->q_angle,*p->or);
	}
}
