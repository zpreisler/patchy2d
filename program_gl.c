#include <stdio.h>
#include <stdlib.h>
#include <utils.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#define MAX_SOURCE_SIZE 8096
char *read_file(char *source){
	char *source_str;
	int source_size;
	FILE *f;
	f=open_file(source,"r");
	if(!f)error("Failed to open file %s",source);
	source_str=(char*)alloc(MAX_SOURCE_SIZE);
	source_size=fread(source_str,1,MAX_SOURCE_SIZE,f);
	if(!source_size)error("Failed to read file %s",source);
	close_file(f);
	return source_str;
}
void print_shader_log(unsigned int obj){
	int n;
	char *log;
	glGetShaderiv(obj,GL_INFO_LOG_LENGTH,&n);
	if(n>0){
		log=(char*)alloc(n);
		glGetShaderInfoLog(obj,n,NULL,log);
		error("\n%s",log);
		free(log);
	}
	return;
}
void print_program_log(unsigned int obj){
	int n;
	char *log;
	glGetProgramiv(obj,GL_INFO_LOG_LENGTH,&n);
	if(n>0){
		log=(char*)alloc(n);
		glGetProgramInfoLog(obj,n,NULL,log);
		error("%s",log);
		free(log);
	}
	return;
}
unsigned int load_compile_shader(char *name,GLenum shadertype){
	const char *source;
	unsigned int shader;
	source=read_file(name);
	shader=glCreateShader(shadertype);
	glShaderSource(shader,1,&source,NULL);
	glCompileShader(shader);
	print_shader_log(shader);
	return shader;
}
unsigned int create_program(char *path_vertex_shader,char *path_geometry_shader,char *path_fragment_shader){
	unsigned int program;
	unsigned int vertex_shader,geometry_shader,fragment_shader;
	//Load shaders
	vertex_shader=load_compile_shader(path_vertex_shader,GL_VERTEX_SHADER);
	if(path_geometry_shader){
		geometry_shader=load_compile_shader(path_geometry_shader,GL_GEOMETRY_SHADER);
	}
	fragment_shader=load_compile_shader(path_fragment_shader,GL_FRAGMENT_SHADER);
	//Create program
	program=glCreateProgram();
	//Attach shaders
	glAttachShader(program,vertex_shader);
	if(path_geometry_shader){
		glAttachShader(program,geometry_shader);
	}
	glAttachShader(program,fragment_shader);
	//Link program
	glLinkProgram(program);
	//glUseProgram(program);
	//Clean up
	glDeleteShader(vertex_shader);
	if(path_geometry_shader){
		glDeleteShader(geometry_shader);
	}
	glDeleteShader(fragment_shader);
	//Return program
	return program;
}
void glsl_info(){
	printf("VENDOR: %s\n",glGetString(GL_VENDOR));
	printf("RENDERER: %s\n",glGetString(GL_RENDERER));
	printf("VERSION: %s\n",glGetString(GL_VERSION));
	printf("GLSL: %s\n",glGetString(GL_SHADING_LANGUAGE_VERSION));
}
