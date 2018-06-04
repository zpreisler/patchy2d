#version 450 core
layout(location=0) in vec2 position;

uniform mat2 proj_matrix;
uniform mat2 view_matrix;

void main(){
	vec2 x=position;
	x=(2.0*(x-view_matrix[0])/(view_matrix[1]-view_matrix[0])-1.0)*proj_matrix;
	gl_Position=vec4(x,0.0,1.0);
}
