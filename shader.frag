#version 450 core
out vec4 out_color;

in vec4 col;
in vec2 uv;

void main(){
	float dist=dot(uv,uv);
	if(dist<0.8){
		out_color=vec4(1.0,0.0,0.0,0.1);
	}
	else if(dist<1.0){
		out_color=vec4(col.xyz,0.9);
	}
	else{
		out_color=vec4(1.0,0.0,1.0,0.0);
	}
}
