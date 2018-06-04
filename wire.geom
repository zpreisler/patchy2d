#version 450
layout (triangles) in;
layout (triangle_strip,max_vertices=3) out;

out vec2 uv;
out vec3 dist;

void main(){
	vec2 p0=gl_in[0].gl_Position.xy;
	vec2 p1=gl_in[1].gl_Position.xy;
	vec2 p2=gl_in[2].gl_Position.xy;

	vec2 v0=p2-p1;
	vec2 v1=p2-p0;
	vec2 v2=p1-p0;
	float area=abs(v1.x*v2.y-v1.y*v2.x);
	
	uv=vec2(-1.0,-1.0);
	dist=vec3(area/length(v0),0,0);
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	uv=vec2(-1.0,1.0);
	dist=vec3(0,area/length(v1),0);
	gl_Position = gl_in[1].gl_Position;
	EmitVertex();

	uv=vec2(1.0,-1.0);
	dist=vec3(0,0,area/length(v2));
	gl_Position = gl_in[2].gl_Position;
	EmitVertex();

	EndPrimitive();
}
