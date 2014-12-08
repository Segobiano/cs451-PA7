#pragma once

#include "GL/gli.h"
#include "model.h"
#include <list>

struct Ray
{
	Point3d  o; //origin
	Vector3d v; //direction 
};

class RayTracer
{
public:
	Vector3d inter;
	double mImage_w;
	double mImage_h;

	GLdouble* modelmax;
	GLdouble* project;
	GLint* view;
	GLdouble* scx;
	GLdouble* scy;
	GLdouble* scz;
	GLdouble* scxf;
	GLdouble* scyf;
	GLdouble* sczf;

	double xm_min;
	double ym_min;

	
	double ab;
	double bc;
	double ca;

	double farside;
	double shortside;
	double middleside;


	std::vector<Vector3d> lightbox;

	float camx;
	float camy;
	float camz;
	//constructor
	RayTracer(list<object3D*>& models);

	//render an image with "nray_per_pixel" rays created per pixel
	void render(unsigned int img_w, unsigned int img_h, unsigned int nray_per_pixel);

	//save rendered image to file
	bool save2file(const string& filename);

private:

	// render an image with "nray_per_pixel" rays created per pixel
	// create this ray randomly in the given pixel (x,y)
	Ray create_a_random_ray(unsigned int x, unsigned int y);
	
	Vector3d randPoint();
	//returns a model and triangle that intersect ray r
	//return pair<NULL,NULL> if no intersection is found
	pair<object3D *, triangle *> intersect(Ray r);

	//returns a triangle in model m that intersect ray r
	//return NULL if no intersection is found
	triangle * intersect(model& m, Ray r);

	//returns a triangle in model m that make closest intersection with ray r
	//return NULL if no intersection is found
	triangle * closest_intersect(model& m, Ray r);

	// determine if there is an intersection between triangle t and ray r
	// return true if there is an intersection and store the intersection in x
	// return false otherwise and x is undefined in this case
	bool intersect(model& m, triangle * t, Ray r, Point3d& x);

	// returns true if the sphere intersects ray r
	// return false if no intersection is found
	// x is the location of intersection
	bool intersect(sphere& s, Ray r, Point3d& x);

	//check if a point p is in shadow
	double inshadow();

	//
	// determine the color of ray r
	// 
	Vector3d raycolor(const Ray& r, int depth);

	//
	// determine the color of ray r, by analizing the intersection between t and r 
	// 
	Vector3d raycolor(model& m, triangle * t, const Ray& r, int depth);

	//
	// determine the color of ray r, by analizing the intersection between t and r 
	// 
	Vector3d raycolor(sphere& s, const Ray& r, int depth);

	list<object3D*>& m_models;
	list<object3D*>  m_lights; //this is a subset of m_models
	vector<object3D*> lightSrc;

	vector< vector< Vector3d > > m_img; //rendered image

public:

	//for debugging only
	list<Ray> all_rays;
	list<Point3d> intersection_points;
};

