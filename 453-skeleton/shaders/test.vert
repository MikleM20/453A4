#version 330 core
layout (location = 0) in vec3 pos;
layout (location = 1) in vec2 texCoord;
layout (location = 2) in vec3 normal;


uniform mat4 M;
uniform mat4 V;
uniform mat4 P;
uniform mat4 transform;
uniform int sunObject;

out vec3 fragPos;
out vec2 tc;
out vec3 n;
flat out int isSun;

void main() {
	isSun  = sunObject;
	fragPos = pos;
	tc = texCoord;
	n = normal;
	gl_Position = P * V *M *transform*vec4(pos, 1.0);
}
