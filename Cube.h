#pragma once
#include <cmath>
#include <vector>
#include "Vec3.h"

class Cube{
private:
	Vec3f m_p1;
	Vec3f m_p2;
public:
	Cube(){};
	Cube(Vec3f p1,Vec3f p2){
		m_p1=p1; //xmin ymax zmin
		m_p2=p2; //xmax ymin zmax

	}
	Vec3f getp1() const {return m_p1;};
	Vec3f getp2() const {return m_p2;};


	bool containsVertex(Vec3f vertex){
		return(vertex[0]>=m_p1[0] && vertex[0]<=m_p2[0] && vertex[1]<=m_p1[1] && vertex[1]>=m_p2[1] && vertex[2]>=m_p1[2] && vertex[2]<=m_p2[2]);
	}


	bool containsTriangle(Vec3f p1,Vec3f p2,Vec3f p3); //deprecated



};