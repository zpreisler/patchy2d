#version 450
layout (points) in;
layout (triangle_strip,max_vertices=4) out;

in vec4 vertex_color[1];
out vec4 geom_color;
out vec2 uv;

uniform mat2 proj_matrix;
uniform mat2 view_matrix;

void main(){
	vec2 a=vec2(1.0,1.0)*proj_matrix/(view_matrix[1]-view_matrix[0]);
	geom_color=vertex_color[0];
	uv=vec2(-1.0,-1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	geom_color=vertex_color[0];
	uv=vec2(-1.0,1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	geom_color=vertex_color[0];
	uv=vec2(1.0,-1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	geom_color=vertex_color[0];
	uv=vec2(1.0,1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	EndPrimitive();
}
