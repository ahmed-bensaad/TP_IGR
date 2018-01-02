#include "Cube.h"



bool Cube::containsTriangle(Vec3f p1,Vec3f p2,Vec3f p3){
	unsigned int i=0;
	std::vector<Vec3f> vertices;
	vertices.push_back(p1);
	vertices.push_back(p2);
	vertices.push_back(p3);
	for(std::vector<Vec3f>::iterator it = vertices.begin(); it != vertices.end(); ++it){
	if (containsVertex(*it)){
		i+=1;
		}
	}
	return(i>=2);


}
