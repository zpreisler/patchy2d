#version 450 core
layout(location=0) in dvec2 position;
layout(location=1) in vec4 color;

uniform mat2 proj_matrix;
uniform mat2 view_matrix;

out vec4 vertex_color;

void main(){
	dvec2 x=position;
	x=x*proj_matrix;
	gl_Position=vec4(x,0.0,1.0);
	vertex_color=color;
}
