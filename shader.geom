#version 450
layout (points) in;
layout (triangle_strip,max_vertices=16) out;

in vec4 vertex_color[];
out vec4 geom_color;
out vec2 uv;

uniform mat2 proj_matrix;
uniform mat2 view_matrix;

void draw_sphere(vec4 position,vec4 color,vec2 a,vec2 b){
	vec4 x=vec4(b,0.0,0.0);
	geom_color=color;
	uv=vec2(-1.0,-1.0);
	gl_Position=position+x+vec4(a*uv,0.0,0.0);
	EmitVertex();

	geom_color=color;
	uv=vec2(-1.0,1.0);
	gl_Position=position+x+vec4(a*uv,0.0,0.0);
	EmitVertex();

	geom_color=color;
	uv=vec2(1.0,-1.0);
	gl_Position=position+x+vec4(a*uv,0.0,0.0);
	EmitVertex();

	geom_color=color;
	uv=vec2(1.0,1.0);
	gl_Position=position+x+vec4(a*uv,0.0,0.0);
	EmitVertex();

	EndPrimitive();
}

void main(){
	vec2 a,b;
	vec4 p=gl_in[0].gl_Position;
	vec4 c=vertex_color[0];

	a=vec2(0.5,0.5)*proj_matrix;

	b=vec2(0.0,0.0);
	draw_sphere(p,c,a,b);
	
	b=vec2(-1.0,0.0);
	draw_sphere(p,c,a,b);

	b=vec2(0.0,-1.0);
	draw_sphere(p,c,a,b);

	b=vec2(-1.0,-1.0);
	draw_sphere(p,c,a,b);
}
