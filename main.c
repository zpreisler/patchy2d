#include <SDL2/SDL.h>
#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES
#include <utils.h>
#define MAX_SOURCE_SIZE (0x100000)
unsigned int vao[10];
unsigned int vbo[10];
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
	geometry_shader=load_compile_shader(path_geometry_shader,GL_GEOMETRY_SHADER);
	fragment_shader=load_compile_shader(path_fragment_shader,GL_FRAGMENT_SHADER);
	//Create program
	program=glCreateProgram();
	//Attach shaders
	glAttachShader(program,vertex_shader);
	glAttachShader(program,geometry_shader);
	glAttachShader(program,fragment_shader);
	//Link program
	glLinkProgram(program);
	glUseProgram(program);
	//Clean up
	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	return program;
}
void initialize(){
	unsigned int program;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float v[]={-0.5,-0.5,0.5,-0.5,0.0,0.5};
	float v2[]={-0.5,0.5,0.5,0.5,0.0,-0.5};
	float color[]={
		0.5,1.0,0.5,1.0,
		0.5,1.0,0.75,1.0,
		0.5,1.0,0.35,1.0
	};
	program=create_program("shader.vert","shader.geom","shader.frag");

	glGenVertexArrays(1,vao);
	glGenBuffers(2,vbo); //Generate buffers

	glBindVertexArray(vao[0]);
	glBindBuffer(GL_ARRAY_BUFFER,vbo[0]); //Select buffer
	glBufferData(GL_ARRAY_BUFFER,sizeof(v),v,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER,vbo[1]); //Select buffer
	glBufferData(GL_ARRAY_BUFFER,sizeof(color),color,GL_STREAM_DRAW); //Write into the buffer
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,0);
	glEnableVertexAttribArray(1);
}
void display(){
	glClearColor(1.0,1.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glBindVertexArray(vao[0]);
	glDrawArrays(GL_POINTS,0,3);
}
void glsl_info(){
	printf("VENDOR: %s\n",glGetString(GL_VENDOR));
	printf("RENDERER: %s\n",glGetString(GL_RENDERER));
	printf("VERSION: %s\n",glGetString(GL_VERSION));
	printf("GLSL: %s\n",glGetString(GL_SHADING_LANGUAGE_VERSION));
}


int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window *window=SDL_CreateWindow("My SDL Empty Window",
			SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED,640,640,
			SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);

	SDL_GLContext glcontext=SDL_GL_CreateContext(window);
	GLenum glew_status = glewInit();

	glsl_info();
	initialize();

	display();
	SDL_GL_SwapWindow(window);
	SDL_Delay(1000);

	SDL_Quit();
	return 0;
}
