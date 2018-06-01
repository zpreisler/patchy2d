#version 450
layout (points) in;
layout (triangle_strip,max_vertices=4) out;

in vec4 vcolor[];
out vec4 col;
out vec2 uv;

void main(){
	vec2 a=vec2(0.1,0.1);

	col=vcolor[0];
	uv=vec2(-1.0,-1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	col=vcolor[0];
	uv=vec2(-1.0,1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	col=vcolor[0];
	uv=vec2(1.0,-1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	col=vcolor[0];
	uv=vec2(1.0,1.0);
	gl_Position = gl_in[0].gl_Position + vec4(a*uv,0.0,0.0);
	EmitVertex();

	EndPrimitive();
}
