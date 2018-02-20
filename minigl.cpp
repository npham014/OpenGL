/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h



/**
 * Global Variables
 */
vec3 color;
MGLpoly_mode currentPrim;
MGLmatrix_mode currentMatrixMode;
mat4 modelViewMatrix;
mat4 projMatrix;
/**
 * Structs n stuff
 */
struct vertex {
	vec4 pos;
	vec3 vertColor;
};

struct triangle {
	vertex a, b, c;
};

vector<vertex> vertices;
vector<triangle> triangles;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

mat4& getCurrentMatrix() { 
	if(currentMatrixMode == MGL_MODELVIEW) {
		return modelViewMatrix;
	}
	else { //Not sure if we need to account for MGL_TEXTURE or MGL_COLOR yet
		return projMatrix;
	}
}

float triArea(vec2 a, vec2 b, vec2 c) { //Returns the area of the triangle made by vertexes a b and c.
	return abs( (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])) / 2.0);
}

void Rasterize_Triangle(const triangle& tri, int width, int height, MGLpixel* data) {
	float fi = ((tri.a.pos[0]/tri.a.pos[3]) + 1) * width;
	fi /= 2.0;
	fi -= 0.5;
	
	float fj = ((tri.a.pos[1]/tri.a.pos[3]) + 1) * height;
	fj /= 2.0;
	fj -= 0.5;

	vec2 pointa = vec2(fi,fj);

	fi = ((tri.b.pos[0]/tri.b.pos[3]) + 1) * width;
	fi /= 2.0;
	fi -= 0.5;
	
	fj = ((tri.b.pos[1]/tri.b.pos[3]) + 1) * height;
	fj /= 2.0;
	fj -= 0.5;

	vec2 pointb = vec2(fi,fj);
	
	fi = ((tri.c.pos[0]/tri.c.pos[3]) + 1) * width;
	fi /= 2.0;
	fi -= 0.5;
	
	fj = ((tri.c.pos[1]/tri.c.pos[3]) + 1) * height;
	fj /= 2.0;
	fj -= 0.5;
	
	vec2 pointc = vec2(fi, fj);

	float totalArea = triArea(pointa, pointb, pointc);

	for(int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			vec2 I = vec2(i,j);
			float alpha = triArea(I, pointb, pointc)/totalArea;	
			float beta = triArea(pointa, I, pointc)/totalArea;
			float gamma = triArea(pointa, pointb, I)/totalArea;
			if(alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0) {
				//Color pixel: don't really know what im doing but yolo
				data[i+j*width] = Make_Pixel(tri.a.vertColor[0] * 255, tri.a.vertColor[0] * 255, tri.a.vertColor[0] * 255); //This will probably change later. Currently using color of 1st pixel 
			}
		}
	}
						
}

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < height; j++) {
			data[i+j*width] = Make_Pixel(color[0] * 255, color[1] * 255, color[2] * 255); //Set all pixels to global color
		}
	}
	
	for(unsigned int k = 0; k < triangles.size(); k++) {
		Rasterize_Triangle(triangles.at(k), width, height, data);
	}	
	
	triangles.clear();
	
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	currentPrim = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(currentPrim == MGL_TRIANGLES) {
		int temp = vertices.size(); //If there aren't an even 3 vertices, don't count the last 2
		temp -= temp % 3;

		for(int i = 0; i < temp; i += 3) { //Loop through the vertices, every 3 creates a triangle.
			triangle t;
			t.a = vertices.at(i);
			t.b = vertices.at(i+1);
			t.c = vertices.at(i+2);
			triangles.push_back(t);
		}
	}

	if(currentPrim == MGL_QUADS) {
		int temp = vertices.size();
		temp -= temp % 4;
		for(int i = 0; i < temp; i +=4) {
			triangle t1;
			triangle t2;
			t1.a = vertices.at(i);
			t1.b = vertices.at(i+1);
			t1.c = vertices.at(i+2);
			
			t2.a = vertices.at(i);
			t2.b = vertices.at(i+2);
			t2.c = vertices.at(i+3);
			
			triangles.push_back(t1);
			triangles.push_back(t2);
		}
	}
	
	vertices.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	vertex v;
	v.pos = vec4(x,y,0,1);
	v.vertColor = color;
	v.pos = projMatrix * v.pos;
	vertices.push_back(v);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	vertex v;
	v.pos = vec4(x,y,z,1);
	v.vertColor = color;
	//v.pos = modelViewMatrix * projMatrix * v.pos;
	v.pos = projMatrix * v.pos; //Change to the above later	
	vertices.push_back(v);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	currentMatrixMode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	//Note: Why do we pass in MGLfloat instead of just passing the matrix???????
	
	mat4 temp;
	temp.make_zero();
	for(int i = 0; i < temp.cols(); i++) {
		for(int j = 0; j < temp.rows(); j++) {
			temp(i,j) = matrix[j*temp.cols() + i];
		}
	}
	 
	mat4& current = getCurrentMatrix();
	current = temp * current;

	return;	
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	float x = (left + right)/(right - left) * -1;
	float y = (top + bottom)/(top - bottom) * -1;
	float z = (far + near)/(far - near) * -1;
	
	float a = 2.0 / (right - left);
	float b = 2.0 / (top - bottom);
	float c = 2.0 / (far - near) * -1;
	
	
	mat4 orthoMatrix = {a, 0.0, 0.0, 0.0,
			    b, 0.0, 0.0, 0.0,
			    c, 0.0, 0.0, 0.0,
			    x, y, z, 1.0};

	mglMultMatrix(orthoMatrix.values);	
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	color = vec3(red, green, blue);
}
