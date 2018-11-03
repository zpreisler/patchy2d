#version 450 core
out vec4 out_color;

in vec4 geom_color;
in vec2 uv;

void main(){
	float dist=dot(uv,uv);
	float delta=fwidth(dist);
	float alpha=1.0-smoothstep(1.0-2*delta,1.0,dist);
	out_color=geom_color*alpha;

	//if(dist<1.0){
	//	out_color=geom_color;
	//}
	//else{
	//	out_color=vec4(1.0,1.0,1.0,0.0);
	//}
}
