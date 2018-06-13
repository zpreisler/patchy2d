#version 450 core
layout(location=0) in vec2 position;
layout(location=1) in vec4 color;

uniform mat2 proj_matrix;
uniform mat2 view_matrix;
uniform float uy;

out vec4 vertex_color;

void main(){
	vec2 x=vec2(position.x+position.y*uy,position.y);
	//x=x*proj_matrix*1.0-vec2(0.5,0.5);
	x=x*proj_matrix;
	gl_Position=vec4(x,0.0,1.0);
	vertex_color=color;
}
