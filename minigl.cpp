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
mat4 modelViewMatrix = {1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1};
mat4 projMatrix = {1, 0, 0, 0,
		   0, 1, 0, 0,
		   0, 0, 1, 0,
                   0, 0, 0, 1};
const mat4 identity = {1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0, 1, 0,
		 0, 0, 0, 1};

vector<mat4> projStack;
vector<mat4> modViewStack;
vector<vector<float>> zBuffer;//zBuffer.at(1).at(2) denotes the pixel at x = 1, y = 2. ie zBuffer is width by height

//mat4* currentTop = &projStack;

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
	if(currentMatrixMode == MGL_PROJECTION) {
		return projMatrix;
	}
	else { //Not sure if we need to account for MGL_TEXTURE or MGL_COLOR yet
		return modelViewMatrix;
	}
}

void setCurrentMatrix(mat4 matrix) {
	getCurrentMatrix() = matrix;
	
}
float triArea(vec2 a, vec2 b, vec2 c) { //Returns the area of the triangle made by vertexes a b and c.
	//return abs( (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])) / 2.0);
	return a[0]*(b[1]-c[1]) + a[1]*(c[0]-b[0]) + (b[0]*c[1] - b[1]*c[0]);
}

void Rasterize_Triangle(const triangle& tri, int width, int height, MGLpixel* data) {
	
	float fi = (tri.a.pos[0]/tri.a.pos[3] + 1) * width;
	fi /= 2.0;
	fi -= 0.5;
	
	float fj = (tri.a.pos[1]/tri.a.pos[3] + 1) * height;
	fj /= 2.0;
	fj -= 0.5;

	vec2 pointa = vec2(fi,fj);

	fi = (tri.b.pos[0]/tri.b.pos[3] + 1) * width;
	fi /= 2.0;
	fi -= 0.5;
	
	fj = (tri.b.pos[1]/tri.b.pos[3] + 1) * height;
	fj /= 2.0;
	fj -= 0.5;

	vec2 pointb = vec2(fi,fj);
	
	fi = (tri.c.pos[0]/tri.c.pos[3] + 1) * width;
	fi /= 2.0;
	fi -= 0.5;
	
	fj = (tri.c.pos[1]/tri.c.pos[3] + 1) * height;
	fj /= 2.0;
	fj -= 0.5;
	
	vec2 pointc = vec2(fi, fj);
	
	/*cout << "a: " << pointa[0] << pointa[1] << endl;
	cout << "b: " << pointb[0] << pointb[1] << endl;
	cout << "c: " << pointc[0] << pointc[1] << endl;
	*/
	float totalArea = triArea(pointa, pointb, pointc);
	//cout << "triangles size:" << triangles.size() << endl;
	for(int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			
			vec2 I = vec2(i,j);
			float alpha = triArea(I, pointb, pointc)/totalArea;	
			float beta = triArea(pointa, I, pointc)/totalArea;
			float gamma = triArea(pointa, pointb, I)/totalArea;
			/*cout << "alpha: " << alpha << endl;
			cout << "beta: " << beta << endl;
			cout << "gamma: " << gamma << endl;*/
			if(alpha >= 0 && beta >= 0 && gamma >= 0) {
				float zDepth = (alpha * tri.a.pos[2]/tri.a.pos[3]) + (beta * tri.b.pos[2]/tri.b.pos[3]) + (gamma * tri.c.pos[2]/tri.c.pos[3]);
				if(zBuffer.at(i).at(j) > zDepth) {
					float k = (alpha / tri.a.pos[3]) + (beta / tri.b.pos[3]) + (gamma / tri.c.pos[3]);
					alpha = (alpha / tri.a.pos[3]) / k;
					beta = (beta / tri.b.pos[3]) / k;
					gamma = (gamma / tri.c.pos[3]) / k;
					vec3 finalColor = (tri.a.vertColor * 255 * alpha+ tri.b.vertColor * 255 * beta+ tri.c.vertColor * 255 * gamma);
					data[i+j*width] = Make_Pixel(finalColor[0], finalColor[1], finalColor[2]);
					zBuffer.at(i).at(j) = zDepth;
				}//Z buffer implementation if 
			}//Barycentric coordinate if
		}//inner loop
	}//outer loop
						
}//function

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
	zBuffer.resize(width);//Initialize Z buffer for an unkown width and height
	for(unsigned l = 0; l < width; l++) {
		zBuffer.at(l) = vector<float>(height, 2);
		
	}
	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < height; j++) {
			data[i+j*width] = Make_Pixel(0, 0, 0); //Set all pixels to black
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
	v.pos = projMatrix * modelViewMatrix * v.pos;
	//v.pos = projStack.back() * modViewStack.back() * v.pos;
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
	v.pos = vec4(x,y,z,1.0);
	
	v.vertColor = color;
	v.pos = projMatrix * modelViewMatrix * v.pos;
	//v.pos = projStack.back() * modViewStack.back() * v.pos;
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
	/*if(projStack.size() == 0) {
		projStack.push_back(identity);
	}
	if(modViewStack.size() == 0) {
		modViewStack.push_back(identity);
	}*/		
	
	if(currentMatrixMode == MGL_PROJECTION) {
		projStack.push_back(getCurrentMatrix());
	}
	else { //Again assuming that MGLTEXTURE and MGLCOLOR aren't needed yet
		modViewStack.push_back(getCurrentMatrix());
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	
	if(currentMatrixMode == MGL_PROJECTION) {
		if(projStack.size() == 0) {
			cout << "Error: Popping Empty Stack" << endl;
			return;
		}
		setCurrentMatrix(projStack.back());
		projStack.pop_back();
	}	
	else {
		if(modViewStack.size() == 0) {
			cout << "Error: Popping Empty Stack." << endl;
			return;
		}
		setCurrentMatrix(modViewStack.back());
		modViewStack.pop_back();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	setCurrentMatrix(identity);
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
	mat4 temp;
	temp.make_zero();
	for(int i = 0; i < temp.cols(); i++) {
		for(int j = 0; j < temp.rows(); j++) {
			temp(i,j) = *(matrix + (i + j * temp.cols()));
		}
	}
	
	setCurrentMatrix(temp);
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
	mat4 temp;
	temp.make_zero();
	for(int i = 0; i < temp.cols(); i++) {
		for(int j = 0; j < temp.rows(); j++) {
			temp(i,j) = *(matrix + (i + j * temp.cols()));
		}
	}
	 
	/*mat4& current = getCurrentMatrix();
	current =  current * temp;*/
	
	getCurrentMatrix() = getCurrentMatrix() * temp;

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
	mat4 transMatrix = {1.0, 0.0, 0.0, 0.0,
			    0.0, 1.0, 0.0, 0.0,
			    0.0, 0.0, 1.0, 0.0,
			    x, y, z, 1.0};

	mglMultMatrix(transMatrix.values);
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
	float c = cos(angle * M_PI / 180.0);
	float s = sin(angle * M_PI / 180.0);
	
	//c *= -1;
	//s *= -1;
	vec3 point = vec3(x,y,z);
	point = point.normalized();
	
	float x1 = point[0];
	float y1 = point[1];
	float z1 = point[2];	
	
	float a1 = ((x1 * x1) * (1 - c)) + c;
	float a2 = ((y1 * x1) * (1 - c)) + (z1 * s);
	float a3 = ((x1 * z1) * (1 - c)) - (y1 * s);
	float b1 = ((x1 * y1) * (1 - c)) - (z1 * s);
	float b2 = ((y1 * y1) * (1 - c)) + c;
	float b3 = ((y1 * z1) * (1 - c)) + (x1 * s);
	float c1 = ((x1 * z1) * (1 - c)) + (y1 * s);
	float c2 = ((y1 * z1) * (1 - c)) - (x1 * s);
	float c3 = ((z1 * z1) * (1 - c)) + c;
	
	mat4 rotMatrix = {a1, a2, a3, 0.0,
			  b1, b2, b3, 0.0,
			  c1, c2, c3, 0.0,
			  0.0, 0.0, 0.0, 1.0};

	mglMultMatrix(rotMatrix.values);	
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat4 scaleMatrix = {x, 0.0, 0.0, 0.0,
			    0.0, y, 0.0, 0.0,
			    0.0, 0.0, z, 0.0,
			    0.0, 0.0, 0.0, 1.0};

	mglMultMatrix(scaleMatrix.values);	
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
	float x = (near * 2.0)/(right - left);
	float y = (near * 2.0)/(top - bottom);
	
	float a = (right + left)/(right - left);
	float b = (top + bottom)/(top - bottom);
	float c = (far + near)/(far - near) * -1.0;
	float d = (2 * far * near)/(far - near) * -1.0;
	
	mat4 frustMatrix = {x, 0.0, 0.0, 0.0,
			    0.0, y, 0.0, 0.0,
			    a, b, c, -1.0,
			    0.0, 0.0, d, 0.0};

	
	mglMultMatrix(frustMatrix.values);		    
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
	float x = (right + left)/(right - left) * -1;
	float y = (top + bottom)/(top - bottom) * -1;
	float z = (far + near)/(far - near) * -1;
	
	float a = 2.0 / (right - left);
	float b = 2.0 / (top - bottom);
	float c = -2.0 / (far - near);
	
	
	mat4 orthoMatrix = {a, 0.0, 0.0, 0.0,
			    0.0, b, 0.0, 0.0,
			    0.0, 0.0, c, 0.0,
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
