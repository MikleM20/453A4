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




struct SceneObject {
	SceneObject(std::string texturePath, GLenum textureInterpolation) :
		texture(texturePath, textureInterpolation),
		position(0.0f, 0.0f, 0.0f),
		scale(1),
		transformationMatrix(1.0f) // This constructor sets it as the identity matrix
	{}

	CPU_Geometry cgeom;
	GPU_Geometry ggeom;
	Texture texture;

	glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);

	// Alternatively, you could represent rotation via a normalized heading vec:
	// glm::vec3 heading;
	glm::vec2 scale; // Or, alternatively, a glm::vec2 scale;
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
void makeSphere() {
	float step = 1.f / (float)(Sections - 1);
	float u = 0.f;

	// Traversing the planes of time and space
	for (int i = 0; i < Sections; i++) {
		float v = 0.f;

		//Traversing the planes of time and space (again)
		for (int j = 0; j < Sections; j++) {
			glm::vec3 vertex = glm::vec3(1 * cos(2.f * M_PI * u) * sin(M_PI * v),
				1 * sin(2.f * M_PI * u) * sin(M_PI * v),
				1 * cos(M_PI * v));

			glm::vec3 normal = glm::vec3(vertex);
			numberOfPoints++;
			points.push_back(vertex);
			normals.push_back(normal);
			textureMap.push_back(glm::vec2(u, v));

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

	SceneObject earth("453-skeleton/textures/2k_earth_daymap.png", GL_LINEAR);
	makeSphere();


	
	CPU_Geometry sphere;
	for (int i = 0; i < numberOfIndexes; i++) {
		int index = indexArray[i];
		//Log::debug("Index");
		earth.cgeom.verts.push_back(points[index]);
		//earth.cgeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
		earth.cgeom.normals.push_back(points[index]);
		earth.cgeom.texCoords.push_back(textureMap[index]);


		sphere.verts.push_back(points[index]);
		//earth.cgeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
		sphere.normals.push_back(points[index]);
		sphere.texCoords.push_back(textureMap[index]);
	}

	updateGPUGeometry(earth.ggeom, earth.cgeom);



	/**
	std::vector<glm::vec3> originQuad;
	originQuad.push_back(glm::vec3{-0.5, 0.5, 0}); // top-left
	originQuad.push_back(glm::vec3{-0.5, -0.5, 0}); // bottom-left
	originQuad.push_back(glm::vec3{0.5, 0.5, 0}); // top-right

	originQuad.push_back(glm::vec3{-0.5, -0.5, 0}); // bottom-left
	originQuad.push_back(glm::vec3{0.5, -0.5, 0}); // bottom-right
	originQuad.push_back(glm::vec3{0.5, 0.5, 0}); // top-right

	

	CPU_Geometry square;
	positiveZFace(originQuad, square);
	positiveXFace(originQuad, square);
	negativeZFace(originQuad, square);
	negativeXFace(originQuad, square);
	positiveYFace(originQuad, square);
	negativeYFace(originQuad, square);

	//square.cols.resize(square.verts.size(), glm::vec3{1.0, 0.0, 0.0});
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);


	for(auto i = square.verts.begin(); i < square.verts.end(); ++i) {
		std::cout << *i << std::endl;
	}

	GPU_Geometry quads;
	updateGPUGeometry(quads, square);
	*/
	Log::debug("Start Loop");
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

		earth.ggeom.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(sphere.verts.size()));

		//quads.bind();
		//glDrawArrays(GL_TRIANGLES, 0, GLsizei(square.verts.size()));

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
