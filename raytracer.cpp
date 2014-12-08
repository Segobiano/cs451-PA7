#include "raytracer.h"
#include <algorithm>
#include <cfloat>
#include "GL/gliLight.h"

//NEW in PA7
#define MAX_RAY_DEPTH 5
#define SHADOW_SAMPLE_SIZE 10

//some helper functions
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int( clamp(x)*255 + .5); }


//constructor
RayTracer::RayTracer(list<object3D*>& models) : m_models(models)
{
	modelmax = new GLdouble[16];
	project = new GLdouble[16];
	view = new GLint[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelmax);
	glGetDoublev(GL_PROJECTION_MATRIX, project);
	glGetIntegerv(GL_VIEWPORT, view);
	scx = new GLdouble(0);
	scy = new GLdouble(0);
	scz = new GLdouble(0);
	scxf = new GLdouble(0);
	scyf = new GLdouble(0);
	sczf = new GLdouble(0);
	camx = gli::getCameraPos()[0];
	camy = gli::getCameraPos()[1];
	camz = gli::getCameraPos()[2];
	//find out where lights are
	m_lights.clear();
	for (list<object3D*>::iterator i = models.begin(); i != models.end(); i++)
	{
		object3D* obj = *i;
		if (obj->mat_emission.normsqr() != 0){
			m_lights.push_back(*i);
			lightSrc.push_back(*i);
		}

	}

	if (m_lights.empty())
	{
		cerr << "! Error: There no light defined. Make sure that some object has has non-zero emission." << endl;
		exit(0);
	}


	triangle t = ((model*)(lightSrc[0]))->tris[0];
	Vector3d aPo = ((model*)(lightSrc[0]))->vertices[t.v[0]].p.get();
	Vector3d bPo = ((model*)(lightSrc[0]))->vertices[t.v[1]].p.get();
	Vector3d cPo = ((model*)(lightSrc[0]))->vertices[t.v[2]].p.get();

	
	std::vector<Vector3d> pointlist;
	pointlist.push_back(aPo);
	pointlist.push_back(bPo);
	pointlist.push_back(cPo);

	xm_min = pointlist[0][0];
	ym_min = pointlist[0][2];
	for (int i = 0; i < pointlist.size(); i++){
		if (pointlist[i][0] < xm_min){
			xm_min = pointlist[i][0];
		}
		if (pointlist[i][1] < ym_min){
			ym_min = pointlist[i][2];
		}
	}
	 ab = (aPo - bPo).norm();
	 bc = (bPo - cPo).norm();
	 ca = (cPo - aPo).norm();

	 farside =max( max(ab, bc), ca);
	 shortside =min( min(ab, bc), ca);
	 middleside = sqrt(farside*farside - shortside*shortside);






}

Vector3d RayTracer::randPoint(){
	//lower the light so it doesnt intersect the ceiling
	//Vector3d output = Vector3d(lightSrc[0]->m_position[0], lightSrc[0]->m_position[1]-10, lightSrc[0]->m_position[2]);
	
	Vector3d output = Vector3d(xm_min + shortside*drand48(), light0_position[1], ym_min + middleside*drand48());
	//std::cout << "xmin: " << xm_min << " zmin: " << ym_min << " shortside: " << shortside*drand48() << " middleside: " << middleside*drand48() << " newX: " << xm_min + middleside*drand48() << " newY: " << ym_min + middleside*drand48() << std::endl;
	/*drand48();
	int negx = rand() % 2;
	int negy = rand() % 2;
	double ranx = drand48();
	double rany = drand48();
	if (negx == 0){
		ranx = -ranx;
	}
	if (negy == 0){
		rany = -rany;
	}
	Vector3d output = Vector3d(lightSrc[0]->m_position[0] + middleside*(ranx / 2), lightSrc[0]->m_position[1] - 5, lightSrc[0]->m_position[2]+middleside*(rany / 2));
	*/
	return output;


}

inline void show_progress_bar(double ratio)
{
	// Show the load bar.
	static const int w = 50;
	int   c = (int)(ratio * w);

	cout << setw(3) << (int)(ratio * 100) << "% [";
	for (int x = 0; x<c; x++) cout << "=";
	for (int x = c; x<w; x++) cout << " ";
	cout << "]\r" << flush;
}

//render an image with "nray_per_pixel" rays created per pixel
void RayTracer::render(unsigned int img_w, unsigned int img_h, unsigned int nray_per_pixel)
{
	mImage_h = img_h;
	mImage_w = img_w;

	//create a black image of size img_w X img_h
	m_img = vector< vector< Vector3d > >(img_h, vector< Vector3d >(img_w, Vector3d(0,0,0) ));
	 
	//initialize the bar
	show_progress_bar(0);

	//generate rays
	for (unsigned int y = 0; y < img_h; y++)
	{

		for (unsigned int x = 0; x < img_w; x++)
		{
			Vector3d color;
			for (unsigned int n = 0; n < nray_per_pixel; n++)
			{
				Ray r=create_a_random_ray(x,y);
				Vector3d rc = raycolor(r, 0);
				color = rc + color;
			}

			m_img[y][x] = color / nray_per_pixel;
		}//end of x

		show_progress_bar(y*1.0 / img_h);

	}//end of y
	show_progress_bar(1.0f);
	cout << endl;
}

// render an image with "nray_per_pixel" rays created per pixel
// create this ray randomly in the given pixel (x,y)
Ray RayTracer::create_a_random_ray(unsigned int x, unsigned int y)
{
	Ray r;

	drand48();
	double xvalue = x / mImage_w*(double)view[2] + (double)drand48();
	double yvalue = (view[3] - y) / mImage_h*(double)view[3] + (double)drand48();

	gluUnProject(xvalue, yvalue, 0, modelmax, project, view, scx, scy, scz);
	gluUnProject(xvalue, yvalue, 1, modelmax, project, view, scxf, scyf, sczf);


	Vector3d pPoint = Vector3d(*scx, *scy, *scz);
	Vector3d pPointf = Vector3d(*scxf, *scyf, *sczf);

	Vector3d camOrg = Vector3d(camx, camy, camz);

	//vector between the camera and the plane
	r.v = ((pPointf - pPoint).normalize());

	//so we start on the front clipping plane
	r.o = pPoint + camOrg;
	//all_rays.push_back(r);

	return r;

}

//returns a model and triangle that intersect ray r
//return pair<NULL,NULL> if no intersection is found
pair<object3D *, triangle *> RayTracer::intersect(Ray r)
{
	double min_dist = FLT_MAX;
	triangle * closest_T = NULL;
	object3D * closest_M = NULL;

	for (list<object3D*>::iterator i = m_models.begin(); i != m_models.end(); i++)
	{
		object3D* obj = *i;

		if (dynamic_cast<model*>(obj)) //this object is a mesh
		{
			model* mesh = dynamic_cast<model*>(obj);
			triangle * t = closest_intersect(*mesh, r);
			if (t != NULL)
			{
				Point3d x;
				intersect(*mesh, t, r, x);
				double dist = (x - r.o).normsqr();
				if (dist < min_dist)
				{
					inter = x.get();
					min_dist = dist;
					closest_T = t;
					closest_M = mesh;
				}
			}
		}//--------------------------------------------------------------------------
		else if (dynamic_cast<sphere*>(obj)) //this object is a sphere
		{
			sphere* ball = dynamic_cast<sphere*>(obj);
			Point3d x;
			if (intersect(*ball, r, x))
			{
				double dist = (x - r.o).normsqr();
				if (dist < min_dist)
				{
					inter = x.get();
					min_dist = dist;
					closest_T = NULL;
					closest_M = ball;
				}
			}
		}//--------------------------------------------------------------------------
	}

	return make_pair(closest_M, closest_T);
}

//
//returns true if the sphere intersects ray r
//return false if no intersection is found
//x is the location of intersection if true is returned.
//
bool RayTracer::intersect(sphere& s, Ray r, Point3d& x)
{

	double a, b, c, discrim,paraT;
	a = r.v[0] * r.v[0] + r.v[1] * r.v[1] + r.v[2] * r.v[2];
	
	b = 2 * r.v[0] * (r.o[0] - s.center[0]) + 2 * r.v[1] * (r.o[1] - s.center[1]) + 2 * r.v[2] * (r.o[2] - s.center[2]);

	c = s.center[0] * s.center[0] + s.center[1] * s.center[1] + s.center[2] * s.center[2] + r.o[0] * r.o[0] + r.o[1] * r.o[1] + r.o[2] * r.o[2] + (-2 * (s.center[0] * r.o[0] + s.center[1] * r.o[1] + s.center[2] * r.o[2]) - s.radius*s.radius);


	//descriminate of the quadratic equation
	discrim = b * b - 4*a*c;


	if (discrim < 0){
		return false;
	}
	else{
		paraT = (-b - sqrt(discrim)) / (2.0 * a);
		if (paraT < 0.0){
			return false;
		}
		x = ((r.o + r.v*paraT));
		return true;
	}

}


//returns a triangle in model m that intersect ray r
//return NULL if no intersection is found
triangle * RayTracer::intersect(model& m, Ray r)
{
	for (unsigned int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if ( intersect(m, &m.tris[i], r, x) )
			return &m.tris[i];
	}

	return NULL;
}

//returns a triangle in model m that make closest intersection with ray r
//return NULL if no intersection is found
triangle * RayTracer::closest_intersect(model& m, Ray r)
{
	double min_dist = FLT_MAX;
	triangle * closest_T=NULL;
	for (unsigned int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if (intersect(m, &m.tris[i], r, x))
		{
			double dist = (x - r.o).normsqr();
			if (dist < min_dist)
			{
				closest_T = &m.tris[i];
				min_dist = dist;
			}
		}//end if
	}//end for i

	return closest_T;
}

// determine if there is an intersection between triangle t and ray r
// return true if there is an intersection and store the intersection in x
// return false otherwise and x is undefined in this case
bool RayTracer::intersect(model& m, triangle * t, Ray r, Point3d& x)
{
	Vector3d v0 = m.vertices[t->v[0]].p.get();
	Vector3d v1 = m.vertices[t->v[1]].p.get();
	Vector3d v2 = m.vertices[t->v[2]].p.get();

	Vector3d mEdg1, mEdg2, mHvec, mSvec, mQvec;
	float mAnum, mFnum, mUnum, mVnum, parat;
	//constructing edges between the vertices
	mEdg1 = v1 - v0;
	mEdg2 = v2 - v0;

	mHvec = r.v%mEdg2;
	mAnum = mEdg1*mHvec;

	//if we are in range of epsilon since 0 could not be hit
	if ((mAnum > -.000001) && (mAnum < .000001)){
		//no intersection
		return false;
	}
	//inverse
	mFnum = 1 / mAnum;

	//both are points so they the - operater is defined for them
	mSvec = r.o - v0;

	//scaled by mFnum
	mUnum = mFnum*(mSvec*mHvec);

	//outside of 0 and 1 then no intersecction
	if (mUnum<0.0 || mUnum>1.0){
		return false;
	}
	//cross product between the edge and the origin and a point on the triangle
	mQvec = mSvec%mEdg1;
	mVnum = mFnum*(r.v*mQvec);

	if (mVnum<0.0 || (mUnum + mVnum)>1.0){
		//it is missing the intersection on the other side
		return false;
	}

	//compute the variable of the parametric
	double paraT = mFnum*(mEdg2*mQvec);

	//intersection happened
	if (paraT > .000001){
		//find the point on of intersection and return true and put the intersection in x

		x = ((r.o + r.v*paraT));
		return true;
	}

	return false;
}


//
// determine the color of ray r
// 
Vector3d RayTracer::raycolor(const Ray& r, int depth)
{
	if (depth >= MAX_RAY_DEPTH) return Vector3d(0,0,0);

	Vector3d color;
	pair<object3D *, triangle *> X = intersect(r);

	//determine the color of this ray 
	if (X.first != NULL && X.second != NULL) //this is mesh intersection
	{
		Vector3d rc = raycolor(*((model*)X.first), X.second, r,depth);
		color = rc + color;
	}
	else if (X.first != NULL) //this is spehere intersection
	{
		Vector3d rc = raycolor(*((sphere*)X.first), r, depth);
		color = rc + color;
	}

	return color;
}

//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(sphere& s, const Ray& r, int depth)
{
	Vector3d color, vv, reflection_color, refraction_color,sNorm,tt;
	Ray nRay;
	double rev;
	//color=s.mat_color;

	sNorm = ((inter - s.center.get()).normalize());
	tt=(lightSrc[0]->m_position - inter).normalize();
	
	color = color + ((s.mat_color ^ lightSrc[0]->mat_emission) * max((sNorm * tt), 0.0));

	//Vector3d temp = inter;
	
	double ins = inshadow();
	if (ins<1){
		return s.mat_color *ins;
	}
	//inter = temp;
	//PA7 TODO: determine the direct color
	//
	//hint: (1) what is the normal direction at the intersection between s and r
	//      (2) what is the color of this sphere? 
	//      (3) where is the light?
	//      (4) is the intersection in shadow? (call inshadow(const Point3d& p))

	//double c1 = -(sNorm*r.v);
	
	if (s.reflectiveness.normsqr() > 0)
	{
		//new ray generated
		
		nRay.v=r.v-(sNorm*((r.v*sNorm)*2.0));
		nRay.o=inter + nRay.v*.001;
		reflection_color = raycolor(nRay, depth+1);
		
		//PA7 TODO: determine the reflection color
		//
		//      (1) generate relection ray 
		//      (2) determine the color of the ray (remember of increase depth by 1)
	}
	
	//inter = temp;
	
	if (s.transparency.normsqr() > 0)
	{

		//entering the shpere
		if ((r.o - s.center).norm() - s.radius < 0.0001){
			rev = s.refractive_index / 1.006;
			//the new direction
			vv = ((sNorm.normalize()  * (rev * (sNorm.normalize() * r.v) - sqrt((1.0f - ((rev * rev) * (1.0f - ((sNorm.normalize() * r.v) * (sNorm.normalize() * r.v)))))))) - (r.v * rev)).normalize();
			
			
			nRay.v = vv;
			nRay.o = inter + nRay.v*.0001;
			//recursive call
			refraction_color = raycolor(nRay, depth+1);
		}
		else{//exiting the sphere
			rev = 1.006 / s.refractive_index;
		
			vv = ((sNorm.normalize()  * (rev * (sNorm.normalize() * r.v) - sqrt((1.0f - ((rev * rev) * (1.0f - ((sNorm.normalize() * r.v) * (sNorm.normalize() * r.v)))))))) - (r.v * rev)).normalize();
			
			//refraction_color = (rev * (sNorm * temp) - ((sqrt(1 - rev * rev * (1.0 - (sNorm * temp) * (sNorm *temp)))))) * sNorm - (rev * temp);
			nRay.v = vv;
			nRay.o = inter + nRay.v*.0001;
			depth++;
			refraction_color = raycolor(nRay, depth+1);
		}
	}
	
		//PA7 TODO: determine the refraction color
		//
		//      (1) generate refraction ray 
		//      (2) determine the color of the ray (remember of increase depth by 1)
	
	//finally gather all the colors
	for (int i = 0; i < 3; i++)
	{
		color[i] = clamp(color[i] + s.reflectiveness[i] * reflection_color[i] + s.transparency[i] * refraction_color[i]);
	}

	return color;
}


//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(model& m, triangle * t, const Ray& r, int depth)
{
	if (m.mat_emission.normsqr() != 0.0){
		return m.mat_emission;
	}
	// The area of a triangle is 
	Vector3d color = m.mat_color;
	//return m.mat_color;
	Vector3d lightPos = Vector3d(light0_position[0], light0_position[1], light0_position[2]);
	//Vector3d lightPos = Vector3d(lightSrc[0]->m_position[0], lightSrc[0]->m_position[1], lightSrc[0]->m_position[2]);
	vertex v0 = m.vertices[t->v[0]];
	vertex v1 = m.vertices[t->v[1]];
	vertex v2 = m.vertices[t->v[2]];
	Vector3d bary;
	//TODO: implement this
	//Interpolate the normal
	float areaABC = t->n * ((Vector3d(v1.p.get()) - Vector3d(v0.p.get())) % (Vector3d(v2.p.get()) - Vector3d(v0.p.get())));
	float areaPBC = t->n * ((Vector3d(v1.p.get()) - inter) % (Vector3d(v2.p.get()) - inter));
	float areaPCA = t->n * ((Vector3d(v2.p.get()) - inter) % (Vector3d(v0.p.get()) - inter));

	bary[0] = areaPBC / areaABC;
	bary[1] = areaPCA / areaABC;
	bary[2] = 1.0f - bary[0] - bary[1];
	//Vector3d weights = getBarycentricCoordinatesAt(Vector3d(pos.get()), Vector3d(v0.p.get()), Vector3d(v1.p.get()), Vector3d(v2.p.get()), t->n);
	Vector3d interNorm = (bary[0] * v0.n) + (bary[1] * v1.n) + (bary[2] * v2.n);
	//Get the diffusion value
	Vector3d fragToLight = (lightPos - inter).normalize();
	float lightval = fragToLight * interNorm;

	//Get the specular value
	Vector3d camOrg = Vector3d(camx, camy, camz);
	Vector3d eyeVec = ((lightPos - inter).normalize() + (camOrg - inter).normalize()).normalize();
	Vector3d specLight = m.mat_specular;
	if (eyeVec*interNorm > 0){
		specLight = specLight* pow((eyeVec*interNorm), m.mat_shininess);
	}
	else{
		specLight = specLight * 0;
	}

	//is the frag facing away from the light?
	if (lightval<0){

		return color * 0;
	}
	//the intersection is in shadow
	double ins = inshadow();
	if (ins<1){
		return m.mat_color *ins;
	}
	else{
		return color*lightval + specLight;
	}

}

//check if a point p is in shadow
//return 0 is p is not in the shadow
//a value between 0 and 1 to indicate the "probability"
//that p is in full darkness
double RayTracer::inshadow()
{
	int sample = 0;
	Vector3d temp = inter;
	for (int i=0; i < SHADOW_SAMPLE_SIZE; i++){
		inter = temp;
		Ray lightRay;
		//random location in the light
		//Vector3d lightPos = Vector3d(lightSrc[0]->m_position[0], lightSrc[0]->m_position[1]-5, lightSrc[0]->m_position[2]);
		Vector3d lightPos = randPoint();
		////////////////////////////////////////////
		lightRay.v = (lightPos - inter).normalize();
		//move of the object
		lightRay.o = inter + lightRay.v*.0001;
		pair<object3D *, triangle *> X = intersect(lightRay);
		
		//intersected with either a sphere or a tri
		if (X.first != NULL){
			/*if (X.first->mat_emission.normsqr() != 0.0){
				sample++;
			}*/
			if ((lightRay.o - inter).norm() < (lightRay.o - lightPos).norm()){
				sample++;
			}
		}
	}
	inter=temp;
	return ((double)SHADOW_SAMPLE_SIZE-sample) / (double)SHADOW_SAMPLE_SIZE;
}

//save rendered image to file
bool RayTracer::save2file(const string& filename)
{
	FILE *f = fopen(filename.c_str(), "w");         // Write image to PPM file.

	int h = m_img.size();
	if (h == 0) return true; //nothing to save...

	int w = m_img.front().size();
	
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
			fprintf(f, "%d %d %d ", toInt(m_img[i][j][0]), toInt(m_img[i][j][1]), toInt(m_img[i][j][2]));

		show_progress_bar(i*1.0 / h);
	}
	show_progress_bar(1.0f);
	fclose(f);
}


