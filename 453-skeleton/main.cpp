#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"
#include "Camera.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

#define Sections 50

int numberOfPoints = 0;
int numberOfIndexes = 0;
std::vector<glm::vec3> points;
std::vector<glm::vec3> normals;
std::vector<glm::vec2> textureMap;
GLuint indexArray[50*50*6];

int numberOfPoints2 = 0;
int numberOfIndexes2 = 0;
std::vector<glm::vec3> points2;
std::vector<glm::vec3> normals2;
std::vector<glm::vec2> textureMap2;
GLuint indexArray2[50 * 50 * 6];
bool play = true;




struct SceneObject {
	SceneObject(std::string texturePath, GLenum textureInterpolation, float tiltIn, float orbitIn, float initialThetaIn, float spinIn, float scaleIn) :
		texture(texturePath, textureInterpolation),
		position(0.0f, 0.0f, 0.0f),
		orbit(orbitIn),
		tilt(tiltIn),
		spin(spinIn),
		scale(scaleIn),
		initialTheta(initialThetaIn), 
		transformationMatrix(1.0f) // This constructor sets it as the identity matrix
	{}

	CPU_Geometry cgeom;
	GPU_Geometry ggeom;
	Texture texture;
	float initialTheta;
	float orbit;
	float spin;
	float tilt;
	glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);

	// Alternatively, you could represent rotation via a normalized heading vec:
	// glm::vec3 heading;
	float scale; // Or, alternatively, a glm::vec2 scale;
	glm::mat4 transformationMatrix;
};


// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	//gpuGeom.setCols(cpuGeom.cols);
	gpuGeom.setTexCoords(cpuGeom.texCoords);
	gpuGeom.setNormals(cpuGeom.normals);
}

// EXAMPLE CALLBACKS
class Assignment4 : public CallbackInterface {

public:
	Assignment4() : camera(0.0, 0.0, 2.0), aspect(1.0f) {
	}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS && key == GLFW_KEY_SPACE) {
			if (play) {
				play = false;
			}
			else if (!play) {
				play = true;
			}
		}
	}
	virtual void mouseButtonCallback(int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_RIGHT) {
			if (action == GLFW_PRESS) {
				rightMouseDown = true;
			} else if (action == GLFW_RELEASE) {
				rightMouseDown = false;
			}
		}
	}
	virtual void cursorPosCallback(double xpos, double ypos) {
		if (rightMouseDown) {
			double dx = xpos - mouseOldX;
			double dy = ypos - mouseOldY;
			camera.incrementTheta(dy);
			camera.incrementPhi(dx);
		}
		mouseOldX = xpos;
		mouseOldY = ypos;
	}
	virtual void scrollCallback(double xoffset, double yoffset) {
		camera.incrementR(yoffset);
	}
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
		aspect = float(width)/float(height);
	}



	void viewPipeline(ShaderProgram &sp) {
		glm::mat4 M = glm::mat4(1.0);
		glm::mat4 V = camera.getView();
		glm::mat4 P = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 1000.f);

		GLint location = glGetUniformLocation(sp, "light");
		glm::vec3 light = camera.getPos();
		//glm::vec3 light =glm::vec3(0, 0, 0);
		glUniform3fv(location, 1, glm::value_ptr(light));

		GLint uniMat = glGetUniformLocation(sp, "M");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(M));
		uniMat = glGetUniformLocation(sp, "V");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(V));
		uniMat = glGetUniformLocation(sp, "P");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(P));

	}

	Camera camera;

private:

	bool rightMouseDown = false;
	float aspect;
	double mouseOldX;
	double mouseOldY;

};
void makeEarth(glm::mat4 &rotateX, glm::mat4 &rotateY) {
	float step = 1.0f / (float)(Sections - 1);
	float u = 0.f;

	// Traversing the planes of time and space
	for (int i = 0; i < Sections; i++) {
		float v = 0.f;
		Log::debug("test");
		//Traversing the planes of time and space (again)
		for (int j = 0; j < Sections; j++) {
			glm::vec3 vertex = glm::vec3(1 * cos(2.f * M_PI * u) * sin(M_PI * v),
				1 * sin(2.f * M_PI * u) * sin(M_PI * v),
				1 * cos(M_PI * v));

			//glm::vec3 normal = glm::vec3(vertex);
			glm::vec3 normal = glm::vec3(vertex.x, vertex.y, vertex.z);
			numberOfPoints++;
			points.push_back(vertex);
			glm::vec4 rotatedNormals = rotateY *rotateX* glm::vec4(normal, 1.0f);
			normals.push_back(glm::vec3(rotatedNormals.x, -rotatedNormals.y, rotatedNormals.z));
			//normals.push_back(glm::vec3(-normal.x, -normal.y, -normal.z));
			float uToPass = 1-u;
			float vToPass = v;
			textureMap.push_back(glm::vec2(uToPass, vToPass));

			v += step;
		}

		u += step;
	}

	for (int i = 0; i < Sections - 1; i++)
	{
		for (int j = 0; j < Sections - 1; j++)
		{
			unsigned int p00 = i * Sections + j;
			unsigned int p01 = i * Sections + j + 1;
			unsigned int p10 = (i + 1) * Sections + j;
			unsigned int p11 = (i + 1) * Sections + j + 1;

			indexArray[numberOfIndexes++] = p00;
			indexArray[numberOfIndexes++] = p10;
			indexArray[numberOfIndexes++] = p01;

			indexArray[numberOfIndexes++] = p01;
			indexArray[numberOfIndexes++] = p10;
			indexArray[numberOfIndexes++] = p11;
		}
	}
}

void makeSun() {
	float step = 1.0f / (float)(Sections - 1);
	float u = 0.f;

	// Traversing the planes of time and space
	for (int i = 0; i < Sections; i++) {
		float v = 0.f;
		Log::debug("test");
		//Traversing the planes of time and space (again)
		for (int j = 0; j < Sections; j++) {
			glm::vec3 vertex = glm::vec3(1 * cos(2.f * M_PI * u) * sin(M_PI * v),
				1 * sin(2.f * M_PI * u) * sin(M_PI * v),
				1 * cos(M_PI * v));

			//glm::vec3 normal = glm::vec3(vertex);
			glm::vec3 normal = glm::vec3(0,0,0);
			numberOfPoints2++;
			points2.push_back(vertex);
			normals2.push_back(glm::vec3(vertex));
			float uToPass = 1 - u;
			float vToPass = v;
			textureMap2.push_back(glm::vec2(uToPass, vToPass));

			v += step;
		}

		u += step;
	}

	for (int i = 0; i < Sections - 1; i++)
	{
		for (int j = 0; j < Sections - 1; j++)
		{
			unsigned int p00 = i * Sections + j;
			unsigned int p01 = i * Sections + j + 1;
			unsigned int p10 = (i + 1) * Sections + j;
			unsigned int p11 = (i + 1) * Sections + j + 1;

			indexArray2[numberOfIndexes2++] = p00;
			indexArray2[numberOfIndexes2++] = p10;
			indexArray2[numberOfIndexes2++] = p01;

			indexArray2[numberOfIndexes2++] = p01;
			indexArray2[numberOfIndexes2++] = p10;
			indexArray2[numberOfIndexes2++] = p11;
		}
	}
}

glm::mat4 makeSpinMat(float angle) {
	return glm::mat4(
		cos(angle*M_PI / 180), 0.0f, -sin(angle * M_PI / 180), 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		sin(angle * M_PI / 180), 0.0f, cos(angle * M_PI / 180), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
		//cos(earth.initialTheta * M_PI / 180), sin(earth.initialTheta * M_PI / 180), 0.0f, 0.0f,
		//-sin(earth.initialTheta * M_PI / 180), cos(earth.initialTheta * M_PI / 180), 0.0f, 0.0f,
		//0.0f, 0.0f, 1.0f, 0.0f,
		//0.0f, 0.0f, 0.0f, 1.0f
	);
}
int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired


	GLDebug::enable();

	// CALLBACKS
	auto a4 = std::make_shared<Assignment4>();
	window.setCallbacks(a4);

	

	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.

	SceneObject earth("453-skeleton/textures/2k_earth_daymap.png", GL_LINEAR, 23, 0, 180, 3,0.05);
	SceneObject sun("453-skeleton/textures/2k_sun.png", GL_LINEAR, 0, 0, 0, 0.5,0.3);
	SceneObject moon("453-skeleton/textures/2k_moon.png", GL_LINEAR, 0, 0, 0, 1,0.01);
	SceneObject bg("453-skeleton/textures/2k_stars.png", GL_LINEAR, 0, 0, 0, 1, 250);
	glm::mat4 base = glm::mat4(1.0f);
	glm::mat4 spinnerY = glm::mat4(
		cos(earth.initialTheta * M_PI / 180), 0.0f, -sin(earth.initialTheta * M_PI / 180), 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		sin(earth.initialTheta * M_PI / 180), 0.0f, cos(earth.initialTheta * M_PI / 180), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
		//cos(earth.initialTheta * M_PI / 180), sin(earth.initialTheta * M_PI / 180), 0.0f, 0.0f,
		//-sin(earth.initialTheta * M_PI / 180), cos(earth.initialTheta * M_PI / 180), 0.0f, 0.0f,
		//0.0f, 0.0f, 1.0f, 0.0f,
		//0.0f, 0.0f, 0.0f, 1.0f
	);
	glm::mat4 spinnerX = glm::mat4{
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f,cos(earth.initialTheta / 2 * M_PI / 180), sin(earth.initialTheta / 2 * M_PI / 180), 0.0f,
		0.0f,-sin(earth.initialTheta / 2 * M_PI / 180), cos(earth.initialTheta / 2 * M_PI / 180), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	
	makeEarth(spinnerY, spinnerX); //spinner X and Y are so that objects in proper orientation

	makeSun();



	
	CPU_Geometry sphere;
	for (int i = 0; i < numberOfIndexes; i++) {
		int index = indexArray[i];
		//Log::debug("Index");
		earth.cgeom.verts.push_back(points[index]);
		//earth.cgeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
		earth.cgeom.normals.push_back(normals[index]);
		earth.cgeom.texCoords.push_back(textureMap[index]);


		sphere.verts.push_back(points[index]);

		sun.cgeom.verts.push_back(points2[index]);
		sun.cgeom.normals.push_back(glm::vec3(-2,-2,-2));	//Doesnt matter since it wont be used
		sun.cgeom.texCoords.push_back(textureMap2[index]);

		moon.cgeom.verts.push_back(points2[index]);
		moon.cgeom.normals.push_back(normals2[index]);
		moon.cgeom.texCoords.push_back(textureMap2[index]);

		bg.cgeom.verts.push_back(points2[index]);
		bg.cgeom.normals.push_back(glm::vec3(-2, -2, -2));	//Doesnt matter since it wont be used
		bg.cgeom.texCoords.push_back(textureMap2[index]);
		
	}

	updateGPUGeometry(earth.ggeom, earth.cgeom);
	updateGPUGeometry(sun.ggeom, sun.cgeom);
	updateGPUGeometry(bg.ggeom, bg.cgeom);

	//earth.transformationMatrix = spinnerY* spinnerX* base;
	glm::mat4 scaleEarth = glm::mat4{
		earth.scale, 0,0,0,
		0,earth.scale,0,0,
		0,0,earth.scale,0,
		0,0,0,1
	};
	glm::mat4 moveEarth = glm::mat4{
		1, 0, 0, 0 ,
		0, 1, 0, 0 ,
		0, 0, 1, 0 ,
		0.75, 0, 0, 1
	};
	glm::mat4 spinEarth = makeSpinMat(earth.spin);

	glm::mat4 tiltEarth = glm::mat4{
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f,cos(earth.tilt / 2 * M_PI / 180), sin(earth.tilt/ 2 * M_PI / 180), 0.0f,
		0.0f,-sin(earth.tilt/ 2 * M_PI / 180), cos(earth.tilt/ 2 * M_PI / 180), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};

	glm::mat4 orbitEarth = makeSpinMat(earth.orbit);
	Log::debug("Start Loop");

	glm::mat4 scaleSun= glm::mat4{
		sun.scale, 0,0,0,
		0,sun.scale,0,0,
		0,0,sun.scale,0,
		0,0,0,1
	};


	glm::mat4 spinSun = makeSpinMat(sun.spin);

	glm::mat4 moveBG = glm::mat4{
		1, 0, 0, 0 ,
		0, 1, 0, 0 ,
		0, 0, 1, 0 ,
		0, -100, 0, 1
	};
	glm::mat4 spinBG = makeSpinMat(bg.spin);
	glm::mat4 scaleBG = glm::mat4{
	bg.scale, 0,0,0,
	0,bg.scale,0,0,
	0,0,bg.scale,0,
	0,0,0,1
	};
	// RENDER LOOP
	while (!window.shouldClose()) {
		Log::debug("Looping");
		glfwPollEvents();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		shader.use();


		a4->viewPipeline(shader);

		if (play) {
			earth.spin += 1;
			earth.orbit += 0.5;
			sun.spin += 0.25;
			bg.spin += 0.01;
		}
		spinEarth = makeSpinMat(earth.spin);
		orbitEarth = makeSpinMat(earth.orbit);
		earth.transformationMatrix = orbitEarth*moveEarth*scaleEarth*spinnerY*spinEarth * tiltEarth* spinnerX *base;
		earth.ggeom.bind();
		GLint myLoc = glGetUniformLocation(shader, "transform");
		glUniformMatrix4fv(myLoc, 1, false, glm::value_ptr(earth.transformationMatrix));
		GLint myLocF = glGetUniformLocation(shader, "sunObject");
		glUniform1i(myLocF, 0);
		earth.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(sphere.verts.size()));
		earth.texture.unbind();


		sun.ggeom.bind();
		spinSun = makeSpinMat(sun.spin);
		sun.transformationMatrix = scaleSun*spinSun*spinnerY*spinnerX*base;
		glUniformMatrix4fv(myLoc, 1, false, glm::value_ptr(sun.transformationMatrix));
		glUniform1i(myLocF, 1);
		sun.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(sphere.verts.size()));
		sun.texture.unbind();

		bg.ggeom.bind();
		spinBG = makeSpinMat(bg.spin);
		bg.transformationMatrix = scaleBG * spinBG  * base;
		glUniformMatrix4fv(myLoc, 1, false, glm::value_ptr(bg.transformationMatrix));
		glUniform1i(myLocF, 1);
		bg.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(sphere.verts.size()));
		bg.texture.unbind();
		
		//quads.bind();
		//glDrawArrays(GL_TRIANGLES, 0, GLsizei(square.verts.size()));

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
