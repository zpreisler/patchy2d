#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <utils.h>
#include "zargs.h"
species *alloc_specie(void){
	int i;
	species *s=(species*)alloc(sizeof(species));
	s->flag=0;
	s->sigma=1.0;
	s->sigma_well=1.0;
	strcpy(s->patch_type,"symmetric");
	s->patch_width=10.0;
	s->npatch=0; //no patches
	s->grand_canonical=0; //disabled
	s->mu=0.0;
	s->interaction_matrix=(double*)alloc(sizeof(double)*64);
	s->angles=(double*)alloc(sizeof(double)*64);
	for(i=0;i<64;i++)*(s->interaction_matrix+i)=1.0;
	for(i=0;i<64;i++)*(s->angles+i)=0.0;
	return s;
}
char *get_subarg(void *value,arg sub,char *argv[],int *i){
	char *err=NULL;
	char *s=*(argv+*i);
	switch(sub.flag){
		case ARG_STRING:
			strcpy((char*)((char*)value+sub.off),s);
			break;
		case ARG_INT:
			*(int*)((char*)value+sub.off)=strtol(s,&err,10);
			if(*err)error("Not an integer '%s'",s);
			break;
		case ARG_DOUBLE:
			*(double*)((char*)value+sub.off)=strtod(s,&err);
			if(*err)error("Not a double '%s'",s);
			break;
	}
	return err;
}
int command_line_subargs(void *value,char *argv[],arg_set *subz,int *i){
	char *msg;
	int j;
	int count;
	arg *sub=subz->args;
	int ns=subz->n;
	do{
		count=0;
		msg=argv[*i+1]; //read following
		if(!msg)return 0;
		else if(*msg++=='-'){
			//Long option
			if(*msg=='-'){
				for(j=0,++msg;j<ns;j++){
					if(sub[j].long_opt&&!strcmp(sub[j].long_opt,msg)){
						(*i)+=2;
						get_subarg(value,sub[j],argv,i);
						count++;
						break;
					}
				}
				if(j==ns)error("Invalid option -- '%s'",msg);
			}
			//Short option
			else{
				for(j=0;j<ns;j++){
					if(sub[j].short_opt==*msg){
						if(*++msg=='\0'){
							(*i)+=2;
							get_subarg(value,sub[j],argv,i);
							count++;
							break;
						}
						else{
							*(argv+(++(*i)))+=2;
							get_subarg(value,sub[j],argv,i);
							count++;
							break;
						}
					}
				}
			}
		}
	}while(count);
	return 0;
}
char *get_arg(arg a,char *argv[],int *i){
	unsigned j;
	int k;
	char *err=NULL;
	char *s=*(argv+*i);
	species *specie=NULL;
	switch(a.flag){
		case ARG_STRING:
			strcpy((char*)a.value,s);
			break;
		case ARG_INT:
			*(int*)a.value=strtol(s,&err,10);
			if(*err)error("Not an integer '%s'",s);
			break;
		case ARG_LONG:
			*(unsigned long long int*)a.value=strtoul(s,&err,10);
			if(*err)error("Not a unsigned long integer '%s'",s);
			break;
		case ARG_DOUBLE:
			*(double*)a.value=strtod(s,&err);
			if(*err)error("Not a double '%s'",s);
			break;
		case ARG_LFVEC:
			for(j=0;j<a.mask-1;(*i)++,j++){
				*((double*)a.value+j)=strtod(argv[*i],&err);
				if(*err)error("Not a double '%s'",argv[*i]);
				if(!argv[*i+1])error("Missing argument '-%c' '--%s'",a.short_opt,a.long_opt);
				if(*argv[*i+1]=='*'){
					for(k=j++;j<a.mask;j++)*((double*)a.value+j)=*((double*)a.value+k);
					break;
				}
			}
			*((double*)a.value+j)=strtod(argv[*i],&err);
			if(*err)error("Not a double '%s'",argv[*i]);
			break;
		case ARG_SPECIES:
			specie=alloc_specie();
			specie->next=*(species**)a.value;
			*(species**)a.value=specie;
			*((int*)((char*)specie+a.off))=strtol(argv[(*i)],&err,10);
			if(*err)error("Not an integer '%s'",argv[*i]);
			command_line_subargs(specie,argv,(arg_set*)a.sub,i);
			break;
		case ARG_UPDATE:
			if(*(species**)a.value){
				specie=*(species**)a.value;
				*(double*)((char*)specie+a.off)=strtod(s,&err);
				if(*err)error("Not a double '%s'",argv[*i]);
				*(species**)a.value=specie->next;
			}
			break;
	}
	return err;
}
int command_line_args(int argc,char *argv[],arg_set argz){
	char *msg;
	int i,j;
	int count=0;
	arg *a=argz.args;
	int n=argz.n;
	for(i=1;i<argc;i++){
		msg=*(argv+i);
		if(*msg++=='-'){
			//Long option
			if(*msg=='-'){
				for(j=0,++msg;j<n;j++){
					if(a[j].long_opt&&!strcmp(a[j].long_opt,msg)){
						i++;
						get_arg(a[j],argv,&i);
						count++;
						break;
					}
				}
				if(j==n)error("Invalid option -- '%s'",msg);
			}
			else if(*msg=='\0')error("Invalid option - ");
			//Short option
			else{
				for(j=0;j<n;j++){
					if(a[j].short_opt==*msg){
						if(*++msg=='\0'){
							if(++i==argc)error("Missing argument '-%c'",a[j].short_opt);
							get_arg(a[j],argv,&i);
							count++;
							break;
						}
						else{
							*(argv+i)+=2;
							get_arg(a[j],argv,&i);
							count++;
							break;
						}
					}
				}
				if(j==n&&*msg)error("Invalid option - '%c'",*msg);
			}
		}
	}
	return count;
}
int skip(FILE *f){
	int c;
	while(isspace(c=getc(f)));
	if(c==COMMENT){
		while((c=getc(f))!='\n');
		return skip(f);
	}
	return ungetc(c,f);
}
int read_string(FILE *f,char *s,int n){
	int c=0,i;
	if(!s||n<0)return -1;
	skip(f);
	for(i=0;i<n&&(c=getc(f))!=EOF&&!isspace(c);i++,s++)*s=c;
	*s='\0';
	if(c)ungetc(c,f);
	return i;
}
int read_stamp(FILE *f,char *s,int n){
	int c=0,i;
	if(!s||n<0)return -1;
	while(isspace(c=getc(f)));
	ungetc(c,f);
	for(i=0;i<n&&(c=getc(f))!=EOF&&!isspace(c);i++,s++)*s=c;
	*s='\0';
	if(c)ungetc(c,f);
	return i;
}
char *get_fsubarg(FILE *f,void *value,arg sub){
	char buffer[NAME_LENGTH];
	char *err=NULL;
	int j,k,nlfvec;
	double *lfvec;
	switch(sub.flag){
		case ARG_STRING:
			read_string(f,(char*)((char*)value+sub.off),NAME_LENGTH);
			break;
		case ARG_INT:
			read_string(f,buffer,NAME_LENGTH);
			*(int*)((char*)value+sub.off)=strtol(buffer,&err,10);
			if(*err)error("Not an integer '%s'",buffer);
			break;
		case ARG_DOUBLE:
			read_string(f,buffer,NAME_LENGTH);
			*(double*)((char*)value+sub.off)=strtod(buffer,&err);
			if(*err)error("Not a double '%s'",buffer);
			break;
		case ARG_LFVEC:
			printf("get");
			nlfvec=*(int*)((char*)value+sub.n);
			lfvec=*(double**)((char*)value+sub.off);
			for(j=0;j<nlfvec;j++){
				read_string(f,buffer,NAME_LENGTH);
				if(j>0&&*buffer=='*'){
					for(k=j-1;j<nlfvec;j++)*(lfvec+j)=*(lfvec+k);
					break;
				}
				*(lfvec+j)=strtod(buffer,&err);
				if(*err)error("Not a double '%s'",buffer);
			}
			break;
		case ARG_LFVEC2:
			nlfvec=*(int*)((char*)value+sub.n);
			lfvec=*(double**)((char*)value+sub.off);
			for(j=0;j<nlfvec*nlfvec;j++){
				read_string(f,buffer,NAME_LENGTH);
				if(j>0&&*buffer=='*'){
					for(k=j-1;j<nlfvec*nlfvec;j++)*(lfvec+j)=*(lfvec+k);
					break;
				}
				*(lfvec+j)=strtod(buffer,&err);
				if(*err)error("Not a double '%s'",buffer);
			}
			break;

	}
	return err;
}
int file_subargs(FILE *f,void *value,arg_set *subz){
	char buffer[NAME_LENGTH];
	int i;
	int count;
	int n=subz->n;
	arg *a=subz->args;
	long int c;
	if(!feof(f)){
		do{
			c=ftell(f);
			count=0;
			if(read_string(f,buffer,NAME_LENGTH)!=-1){
				for(i=0;i<n;i++){
					if(a[i].key){
						if(!strcmp(buffer,a[i].key)){
							get_fsubarg(f,value,a[i]);
							count++;
						}
					}
				}
			}
		}while(count&&!feof(f));
		fseek(f,c,SEEK_SET);
	}	
	return 0;
}

char *get_farg(FILE *f,arg a){
	char buffer[NAME_LENGTH];
	char *err=NULL;
	unsigned j;
	int k;
	species *specie;
	switch(a.flag){
		case ARG_STRING:
			read_string(f,a.value,a.mask);
			break;
		case ARG_INT:
			read_string(f,buffer,NAME_LENGTH);
			*(int*)a.value=strtol(buffer,&err,10);
			if(*err)error("Not an integer '%s'",buffer);
			break;
		case ARG_LONG:
			read_string(f,buffer,NAME_LENGTH);
			*(unsigned long long int*)a.value=strtoul(buffer,&err,10);
			if(*err)error("Not a long '%s'",buffer);
			break;
		case ARG_DOUBLE:
			read_string(f,buffer,NAME_LENGTH);
			*(double*)a.value=strtod(buffer,&err);
			if(*err)error("Not a double '%s'",buffer);
			break;
		case ARG_LFVEC:
			for(j=0;j<a.mask;j++){
				read_string(f,buffer,NAME_LENGTH);
				if(j>0&&*buffer=='*'){
					for(k=j-1;j<a.mask;j++)*((double*)a.value+j)=*((double*)a.value+k);
					break;
				}
				*((double*)a.value+j)=strtod(buffer,&err);
				if(*err)error("Not a double '%s'",buffer);
			}
			break;
		case ARG_SPECIES:
			specie=alloc_specie();
			specie->next=*(species**)a.value;
			*(species**)a.value=specie;
			read_string(f,buffer,NAME_LENGTH);
			*((int*)((char*)specie+a.off))=strtol(buffer,&err,10);
			if(*err)error("Not an integer '%s'",buffer);
			file_subargs(f,specie,(arg_set*)a.sub);
			break;
	}
	return err;
}
int file_args(FILE *f,arg_set argz){
	char buffer[NAME_LENGTH];
	int i;
	int n=argz.n;
	arg *a=argz.args;
	while(!feof(f)){
		if(read_string(f,buffer,NAME_LENGTH)!=-1){
			for(i=0;i<n;i++){
				if(a[i].key){
					if(!strcmp(buffer,a[i].key)){
						get_farg(f,a[i]);
					}
				}
			}
		}
	}	
	return 0;
}
void dump_subargs(FILE *f,void *value,arg_set *subz){
	int i,j;
	int n=subz->n;
	int nlfvec;
	double *lfvec;
	arg *sub=subz->args;
	for(i=0;i<n;i++){
		switch(sub[i].flag){
			case ARG_STRING:
				fprintf(f,"\t%s %s\n",sub[i].key,(char*)((char*)value+sub[i].off));
				break;
			case ARG_INT:
				fprintf(f,"\t%s %d\n",sub[i].key,*(int*)((char*)value+sub[i].off));
				break;
			case ARG_LONG:
				fprintf(f,"\t%s %Lu\n",sub[i].key,*(unsigned long long int*)((char*)value+sub[i].off));
				break;
			case ARG_DOUBLE:
				fprintf(f,"\t%s %.16lf\n",sub[i].key,*(double*)((char*)value+sub[i].off));
				break;
			case ARG_LFVEC:
				nlfvec=*(int*)((char*)value+sub[i].n);
				lfvec=*(double**)((char*)value+sub[i].off);
				fprintf(f,"\t%s ",sub[i].key);
				for(j=0;j<nlfvec;j++)fprintf(f," %lf",*(lfvec+j));
				fputc('\n',f);
				break;
			case ARG_LFVEC2:
				nlfvec=*(int*)((char*)value+sub[i].n);
				lfvec=*(double**)((char*)value+sub[i].off);
				fprintf(f,"\t%s ",sub[i].key);
				for(j=0;j<nlfvec*nlfvec;j++)fprintf(f," %lf",*(lfvec+j));
				fputc('\n',f);
				break;
		}
	}
}
void dump_args(FILE *f,arg_set argz){
	int i;
	unsigned j;
	arg *a=argz.args;
	int n=argz.n;
	species *s;
	for(i=0;i<n;i++){
		switch(a[i].flag){
			case ARG_STRING:
				fprintf(f,"%s %s\n",a[i].key,(char*)a[i].value);
				break;
			case ARG_INT:
				fprintf(f,"%s %d\n",a[i].key,*(int*)a[i].value);
				break;
			case ARG_LONG:
				fprintf(f,"%s %Lu\n",a[i].key,*(unsigned long long int*)a[i].value);
				break;
			case ARG_DOUBLE:
				fprintf(f,"%s %.16lf\n",a[i].key,*(double*)a[i].value);
				break;
			case ARG_LFVEC:
				fprintf(f,"%s",a[i].key);
				for(j=0;j<a[i].mask;j++)fprintf(f," %lf",*((double*)a[i].value+j));
				fputc('\n',f);
				break;
			case ARG_SPECIES:
				s=*((species**)a[i].value);
				while(s){
					fprintf(f,"%s %d\n",a[i].key,*(int*)((char*)s+a[i].off));
					dump_subargs(f,s,a[i].sub);
					s=s->next;
				}
				break;
			default:
				break;
		}
	}
}
void set_args(header *t){
	int i,j;
	int nargs;
	arg *a,*sub;
	arg subargs[]={
		{.short_opt='d',.long_opt="sigma",.key="particle_diameter:",.flag=ARG_DOUBLE,.off=offsetof(species,sigma),.help="Particle diameter"},
		{.short_opt='l',.long_opt="sigma_well",.key="interaction_length:",.flag=ARG_DOUBLE,.off=offsetof(species,sigma_well),.help="Interaction length"},
		{.long_opt="npatch",.key="number_of_patches:",.flag=ARG_INT,.off=offsetof(species,npatch),.help="Number of patches"},
		{.long_opt="interaction_matrix",.key="interaction_matrix:",.flag=ARG_LFVEC2,.off=offsetof(species,interaction_matrix),.n=offsetof(species,npatch),.help="Interaction matrix"},
		{.long_opt="type",.key="patch_type:",.flag=ARG_STRING,.off=offsetof(species,patch_type),.help="Patch type"},
		{.short_opt='w',.long_opt="patch_width",.key="patch_width:",.flag=ARG_DOUBLE,.off=offsetof(species,patch_width),.help="Patch width"},
		{.short_opt='a',.long_opt="patch_angle",.key="patch_angle:",.flag=ARG_DOUBLE,.off=offsetof(species,patch_angle),.help="Patch angle"},
		{.long_opt="angles",.key="angles:",.flag=ARG_LFVEC,.off=offsetof(species,angles),.n=offsetof(species,npatch),.help="Patch angles"},
		{.long_opt="grand_canonical",.key="grand_canonical:",.flag=ARG_INT,.off=offsetof(species,grand_canonical),.help="Grand canonical -- enabled/disabled"},
		{.short_opt='u',.long_opt="mu",.key="chemical_potential:",.flag=ARG_DOUBLE,.off=offsetof(species,mu),.help="Chemical potential"}
	};
	arg_set s={(sizeof(subargs)/sizeof(arg)),subargs},*tmp;
	arg args[]={
		{.short_opt='n',.long_opt="name",.key="name:",.value=t->name,.flag=ARG_STRING,.help="Name",.mask=NAME_LENGTH},
		{.short_opt='s',.long_opt="step",.key="step:",.value=&(t->step),.flag=ARG_LONG,.help="Number of MC steps"},
		{.short_opt='N',.long_opt="specie",.key="specie:",.value=&(t->specie),.flag=ARG_SPECIES,.off=offsetof(species,N),.help="Species",.sub=&s},
		{.long_opt="new_sigma",.value=&(t->update),.flag=ARG_UPDATE,.off=offsetof(species,sigma),.help="Update sigma"}, //update
		{.long_opt="new_width",.value=&(t->update),.flag=ARG_UPDATE,.off=offsetof(species,patch_width),.help="Update width"}, //update
		{.long_opt="new_mu",.value=&(t->update),.flag=ARG_UPDATE,.off=offsetof(species,mu),.help="Update mu"}, //update
		{.long_opt="new_angle",.value=&(t->update),.flag=ARG_UPDATE,.off=offsetof(species,patch_angle),.help="Update angles"}, //update
		//{.long_opt="new_type",.value=&(t->update),.flag=ARG_UPDATE,.off=offsetof(species,patch_type),.help="Update type"}, //update
		{.short_opt='b',.long_opt="box",.key="box:",.value=&(t->box),.flag=ARG_LFVEC,.mask=2,.help="Simulation box lengths x and y"},
		{.short_opt='c',.long_opt="copy",.key="copy:",.value=&(t->copy),.flag=ARG_LFVEC,.mask=2,.help="Copy the simulation box in x and y direction"},
		{.short_opt='e',.long_opt="epsilon",.key="epsilon:",.value=&(t->epsilon),.flag=ARG_DOUBLE,.help="epsilon -- interaction strength"},
		{.long_opt="end_epsilon",.key="end_epsilon:",.value=&(t->end_epsilon),.flag=ARG_DOUBLE,.help="end epsilon -- final interaction strength"},
		{.short_opt='p',.long_opt="pressure",.key="pressure:",.value=&(t->pressure),.flag=ARG_DOUBLE,.help="pressure"},
		{.short_opt='l',.long_opt="lambda",.key="lambda:",.value=&(t->lambda),.flag=ARG_DOUBLE,.help="lambda"},
		{.short_opt='L',.long_opt="lambda_coupling",.key="lambda_coupling:",.value=&(t->lambda_coupling),.flag=ARG_DOUBLE,.help="lambda coupling"},
		{.long_opt="uy",.key="uy:",.value=&(t->uy),.flag=ARG_DOUBLE,.help="uy -- box tilt"},
		{.long_opt="max_move",.key="maximum_displacement:",.value=&(t->max_displacement[0]),.flag=ARG_DOUBLE,.help="Maximum displacement"},
		{.long_opt="max_rot",.key="maximum_rotation:",.value=&(t->max_rotation),.flag=ARG_DOUBLE,.help="Maximum rotation"},
		{.long_opt="max_vol",.key="maximum_volume:",.value=&(t->max_vol),.flag=ARG_DOUBLE,.help="Maximum volume change"},
		{.long_opt="max_uy",.key="maximum_shape:",.value=&(t->max_uy),.flag=ARG_DOUBLE,.help="Maximum shape change"},
		{.long_opt="max_dsigma",.key="maximum_dsigma:",.value=&(t->max_dsigma),.flag=ARG_DOUBLE,.help="Maximum dsigma"},
		{.short_opt='m',.long_opt="mod",.key="mod:",.value=&(t->mod),.flag=ARG_INT,.help="Write modulus"},
		{.long_opt="pmod",.key="pmod:",.value=&(t->pmod),.flag=ARG_INT,.help="Print modulus"},
		{.short_opt='o',.long_opt="optimize",.key="optimize:",.value=&(t->optimize),.flag=ARG_INT,.help="Step size optimization on/off"},
		{.short_opt='v',.long_opt="verbose",.key="verbose:",.value=&(t->verbose),.flag=ARG_INT,.help="Verbose mode on/off"},
		{.long_opt="snapshot",.key="snapshot:",.value=&(t->snapshot),.flag=ARG_INT,.help="Snapshot mode on/off"},
		{.long_opt="seed",.key="seed:",.value=&(t->seed),.flag=ARG_INT,.help="Seed"}
	};
	nargs=(sizeof(args)/sizeof(arg));
	t->argz.n=(sizeof(args)/sizeof(arg));
	t->argz.args=(arg*)alloc(sizeof(arg)*nargs);
	for(i=0;i<nargs;i++){
		a=t->argz.args+i;
		*a=*(args+i);
		if(a->sub){
			tmp=a->sub;
			a->sub=(arg_set*)alloc(sizeof(arg_set));
			((arg_set*)a->sub)->n=tmp->n;
			((arg_set*)a->sub)->args=(arg*)alloc(sizeof(arg)*((arg_set*)a->sub)->n);
			for(j=0;j<((arg_set*)a->sub)->n;j++){
				sub=((arg_set*)a->sub)->args+j;
				*sub=*(tmp->args+j);
			}
		}
	}
}
void print_args(FILE *f,arg_set argz){
	int i;
	arg *a=argz.args;
	int n=argz.n;
	for(i=0;i<n;i++){
		if(a[i].short_opt){
			fprintf(f,"[-%c] [--%s] {%s} %s\n",a[i].short_opt,a[i].long_opt,a[i].key,a[i].help);
		}
		else{
			fprintf(f,"[--%s] {%s} %s\n",a[i].long_opt,a[i].key,a[i].help);
		}
	}
}
void usage(FILE *f,arg_set argz,char *prog){
	fprintf(f,"Usage: %s\n",prog);
	print_args(f,argz);
}
