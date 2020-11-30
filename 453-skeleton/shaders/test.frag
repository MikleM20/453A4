#version 330 core

in vec3 fragPos;
in vec2 tc;
in vec3 n;
flat in int isSun;

uniform vec3 light;
uniform sampler2D sampler;

out vec4 color;

void main() {
	vec4 d = texture(sampler, tc);
	if(isSun < 0.5){
	vec3 test = vec3(-2,-2,-2);
	vec3 test2 = normalize(test);

	vec3 lightDir = normalize(light - fragPos);
	vec3 normal = normalize(n);
    float diff = max(dot(lightDir, normal), 0.0);
	color = vec4(diff * d);
	}
	if(isSun>0.5)
		color = d;
}
