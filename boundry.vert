#version 450 core
layout(location=0) in vec2 position;

uniform mat2 proj_matrix;
uniform mat2 view_matrix;

void main(){
	vec2 x=position;
	x=x*proj_matrix;
	gl_Position=vec4(x,0.0,1.0);
}
