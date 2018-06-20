#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <utils.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "dSFMT.h"
#include "zargs.h"
#include "mm_math.h"
#include "hash.h"
#include "energy.h"
#include "init.h"
#include "patches.h"
#include "alloc.h"
extern dsfmt_t dsfmt;
void pre_set_particle(compound_particle *c,header *t){
	int i;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		*(p)->q=*(c)->q+*(c)->or*i*0.5;
		*(p)->or=*(c)->or;
		boundary(p->q,t->box);
	}
}
void set_particle(compound_particle *c,header *t){
	int i;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		set_patches(p);
		hash_insert(p,t->h1,t->table);
	}
}
int init_configuration_random(header *t){
	unsigned i,k;
	double a;
	compound_particle *c;
	species *s=t->specie;
	while(s){
		for(i=0;i<s->ncompound;i++){
			k=0;
			c=(compound_particle*)(s->c+i);
			a=2.0*M_PI*dsfmt_genrand_open_open(&dsfmt);
			do{
				*(c)->q=t->box*rnd11();
				*(c)->or=sincosa(a);
				pre_set_particle(c,t);
				k++;
				if(k>4096){
					return 1;
				}
			}while(compound_overlap(c,t));
			set_particle(c,t);
		}
		printf(GREEN"-->"CYAN">"BLUE">"RESET" Initialized %d particles with diameter %.3lf"RESET"\n",s->nparticle,s->sigma);
		s=s->next;
	}
	return 0;
}
/*
int init_configuration_random(header *t){
	unsigned i,k;
	double a;
	particle *q;
	species *s=t->specie;
	while(s){
		for(i=0;i<s->nparticle;i++){
			k=0;
			q=(particle*)(s->p+i);
			a=2.0*M_PI*dsfmt_genrand_open_open(&dsfmt);
			do{
				*(q)->q=t->box*rnd11();
				*(q)->or=sincosa(a);
				k++;
				if(k>4096){
					return 1;
				}
			}while(overlap(q,t));
			set_patches(q);
			hash_insert(q,t->h1,t->table);
		}
		printf(GREEN"-->"CYAN">"BLUE">"RESET" Initialized %d particles with diameter %.3lf"RESET"\n",s->nparticle,s->sigma);
		s=s->next;
	}
	return 0;
}
*/
void rev_list(species **head){
	species *prev=NULL;
	species *current=*head;
	species *next;
	while(current){
		next=current->next;
		current->next=prev;
		prev=current;
		current=next;
	}
	*head=prev;
}
void stamp_configuration(FILE *f,time_t *init_time,long long int passed){
	time_t now=time(0);
	char buffer[NAME_LENGTH];
	fprintf(f,CONF_STAMP);
	strftime(buffer,NAME_LENGTH,"%Y-%m-%d %H:%M:%S",localtime(&now));
	fprintf(f," created %s ",buffer);
	strftime(buffer,NAME_LENGTH,"%Y-%m-%d %H:%M:%S",localtime(init_time));
	fprintf(f," started %s run for %.1lfs steps %Ld\n",buffer,difftime(now,*init_time),passed);
}
int save_configuration(char *file,header *t){
	unsigned int i;
	FILE *fconf=open_file2(file,".conf","w");
	species *s;
	particle *p;
	rev_list(&t->specie);
	stamp_configuration(fconf,&t->init_time,t->steps_passed);
	dump_args(fconf,t->argz);
	rev_list(&t->specie);
	s=t->specie;
	while(s){
		for(i=0;i<s->nparticle;i++){
			p=s->p+i;
			fprintf(fconf,"%.12lf %.12lf %.12lf %.12lf\n",(*(p)->q)[0],(*(p)->q)[1],(*(p)->or)[0],(*(p)->or)[1]);
		}
		s=s->next;
	}
	close_file(fconf);
	return 0;
}
int skip_alpha(FILE *f){
	int c;
	while(isspace(c=getc(f)));
	if(isalpha(c)||c==COMMENT){
		while((c=getc(f))!='\n');
		return skip_alpha(f);
	}
	return ungetc(c,f);
}
int read_cline(FILE *f,particle *p,header *t){
	int c __attribute__ ((unused));
	double *q=(double*)p->q;
	double *or=(double*)p->or;
	c=fscanf(f,"%lf %lf %lf %lf\n",q,q+1,or,or+1);
	*(p)->q=*p->q;
	*(p)->or=normalize(*p->or);
	*(p)->q_track=*(p)->q;
	*(p)->q_well=*(p)->q;
	*(p)->or_well=*(p)->or;
	if(overlap(p,t)){
		error("Failed to read a particle position");
	}
	else{
		hash_insert(p,t->h1,t->table);
	}
	return 0;
}
int load_configuration(char *file,header *t){
	char c;
	unsigned i;
	int n;
	long int off;
	species *s=t->specie;
	particle *q;
	FILE *f=open_file(file,"r");
	skip_alpha(f);
	off=ftell(f);
	while((c=fgetc(f))!=EOF){
		if(c=='\n')n++;
	};
	fseek(f,off,SEEK_SET);
	fprintf(stdout,BLUE"Reading configurational file: "RESET"%s\n",file);
	while(s){
		for(i=0;i<s->N;i++){
			q=s->p+i;
			read_cline(f,q,t);
		}
		s=s->next;
	}
	close_file(f);
	return 0;
}
int load_configuration_file(FILE *f,header *t){
	char c;
	unsigned i;
	int n;
	long int off;
	species *s=t->specie;
	particle *q;
	rewind(f);
	skip_alpha(f);
	off=ftell(f);
	while((c=fgetc(f))!=EOF){
		if(c=='\n')n++;
	};
	fseek(f,off,SEEK_SET);
	while(s){
		for(i=0;i<s->N;i++){
			q=s->p+i;
			read_cline(f,q,t);
			set_patches(q);
		}
		s=s->next;
	}
	return 0;
}
void copy_particle(particle *p,particle *q,header *t,__m128d box_copy,__m128d box){
	*q->q=*p->q+box_copy*box;
	*q->or=*p->or;
	set_patches(q);
	hash_insert(q,t->h1,t->table);
}
void copy_configuration(header *t,__m128d box){
	species *s=t->specie;
	unsigned i;
	int x,y;
	int n;
	__m128d box_copy={0.0,0.0};
	particle *p,*q;
	while(s){
		n=s->N;
		for(x=1;x<t->copy[0];x++){
			box_copy[0]=(double)x;
			box_copy[1]=0.0;
			for(i=0;i<s->N;i++){
				p=s->p+i;
				q=s->p+n++;
				copy_particle(p,q,t,box_copy,box);
			}
		}
		s->N=n;
		for(y=1;y<t->copy[1];y++){
			box_copy[0]=0.0;
			box_copy[1]=(double)y;
			for(i=0;i<s->N;i++){
				p=s->p+i;
				q=s->p+n++;
				copy_particle(p,q,t,box_copy,box);
			}
		}
		s->N=n;
		s=s->next;
	}
	t->copy=_mm_set1_pd(1.0);
}
int init_configuration(char *file,header *t){
	FILE *f;
	if(!(f=open_file(file,"r"))){
		fprintf(stdout,RED">"YELLOW">"GREEN">"RESET" Initializing at random "GREEN"<"YELLOW"<"RED"<"RESET"\n");
		if(init_configuration_random(t)){
			error("Random initialization failed, the system is probably too dense");	
		}
		return 0;
	}
	return 0;
}
input_files *find_configurational_files(int argc,char *argv[]){
	char buffer[NAME_LENGTH];
	int i,k;
	int state[argc-1];
	FILE *f;
	input_files *input=(input_files*)alloc(sizeof(input_files));
	for(i=0,input->n=0;i<argc-1;i++){
		state[i]=0;
		if((f=fopen(argv[i+1],"r"))){
			read_stamp(f,buffer,NAME_LENGTH);
			if(!strcmp(buffer,CONF_STAMP)){
				state[i]=1;
				input->n++;
			}
			close_file(f);
		}
	}
	input->file=(char**)alloc(sizeof(char*)*input->n);
	input->f=(FILE**)alloc(sizeof(FILE*)*input->n);
	for(i=0,k=0;i<argc-1;i++){
		if(state[i]){
			input->file[k]=(char*)alloc(sizeof(char)*(strlen(argv[i+1])+1));
			strcpy(input->file[k++],argv[i+1]);
		}
	}
	return input;
}
void print_input_configurational_files(FILE *f,input_files *input){
	int i;
	fprintf(f,CYAN">>>"RESET"Found %d input configurational files\n"RESET,input->n);
	if(input->n){
		for(i=0;i<input->n;i++){
			printf("%s ",input->file[i]);
		}
		fputc('\n',f);
	}
}
int read_input(int argc,char *argv[],input_files *input,header *t){
	FILE *f;
	__m128d box;
	print_input_configurational_files(stdout,input);
	if(!input->n){

		printf(">>>Read command line arguments\n");
		command_line_args(argc,argv,t->argz);
		box=t->box;
		t->box*=t->copy;

		printf(">>>Allocating particles\n");
		alloc_particles(t);
		init_patches(t);
		set_hash(t);

		printf(">>>Initializing random configuration\n");
		if(init_configuration_random(t)){
			error("Failed to initialize configuration");
		}
		copy_configuration(t,box);
	}
	else{
		f=open_file(*input->file,"r");

		file_args(f,t->argz);
		t->update=t->specie;
		command_line_args(argc,argv,t->argz);
		box=t->box;
		t->box*=t->copy;

		alloc_particles(t);
		init_patches(t);
		set_hash(t);

		load_configuration_file(f,t);
		close_file(f);
		copy_configuration(t,box);
	}
	t->max_displacement=_mm_set1_pd(t->max_displacement[0]);
	return 0;
}

