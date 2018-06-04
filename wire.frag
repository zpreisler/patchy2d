#version 450 core
out vec4 out_color;

in vec2 uv;
in vec3 dist;

void main(){
	float d=min(min(dist.x,dist.y),dist.z);
	out_color=vec4(0.0,0.0,1.0,d);
}
