#pragma once
#include <cmath>
#include <vector>
#include "Vec3.h"
#include "Cube.h"
class Mesh;

class Grid{
private:
	Cube m_including;
	unsigned int m_resolution;
	std::vector<Cube> m_cubes;
	std::vector<std::vector<int>> m_contained;
	std::vector<Vec3f> m_representatives;

public:
	Grid(){};
	Grid(Cube including , unsigned int resolution){
		m_including=including;
		Vec3f p1=including.getp1();
		Vec3f p2=including.getp2();
		float xmin=p1[0];
		float xmax=p2[0];
		float ymin=p2[1];
		float ymax=p1[1];
		float zmin=p1[2];
		float zmax=p2[2];
		for (unsigned int x=0;x<resolution;x++){
			for (unsigned int y=0; y<resolution;y++){
				for (unsigned int z=0;z<resolution;z++){
					float x1= xmin+((x+0.f)/resolution)*(xmax-xmin);
					float y1=ymax+((y+0.f)/resolution)*(ymin-ymax);
					float z1=zmin+((z+1.f)/resolution)*(zmax-zmin);
					float x2=xmin+((x+1.f)/resolution)*(xmax-xmin);
					float y2=y1=ymax+((y+1.f)/resolution)*(ymin-ymax);
					float z2=zmin+((z+0.f)/resolution)*(zmax-zmin);
					Vec3f cell_p1=Vec3f(x1,y1,z1);
					Vec3f cell_p2=Vec3f(x2,y2,z2);
					Cube cell=Cube(cell_p1,cell_p2);
					m_cubes.push_back(cell);
				}
			}
		}
	};

	std::vector<Cube> getCubes()const{return m_cubes;};
	std::vector<std::vector<int>> getContained()const{return m_contained;};	
	std::vector<Vec3f> getRep() const{return m_representatives;};

	void ComputeContained(Mesh mesh);	
	void computeRep(Mesh mesh);

};