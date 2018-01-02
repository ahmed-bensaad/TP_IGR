#include "Grid.h"
#include "Mesh.h"


void Grid::ComputeContained(Mesh mesh){
	m_contained.resize(m_cubes.size());
	for(unsigned int i=0; i<mesh.positions().size();i++){
		Vec3f vertex=mesh.positions()[i];
		for (unsigned int j=0;j<m_cubes.size();j++){
			Cube gridCube=m_cubes[j];
			if (gridCube.containsVertex(vertex)){
				m_contained[j].push_back(i);
			}
		}
		}
}

/*void Grid::computeRep(Mesh mesh){
	ComputeContained(mesh);
	m_representatives.resize(m_cubes.size());
	for(std::vector<std::vector<int>>::iterator it = m_contained.begin(); it != m_contained.end(); ++it){
		unsigned int size=(*it).size();
		Vec3f rep=Vec3f();
		for(std::vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2){
			rep+=mesh.positions()[*it2];
		}
		rep=rep/size;
		m_representatives.push_back(rep);
	}
}*/

void Grid::computeRep(Mesh mesh){
	ComputeContained(mesh);
	m_representatives.clear();
	for(unsigned int i=0;i<m_contained.size();i++){
		Vec3f rep=Vec3f(0.0f,0.0f,0.0f);
		unsigned int size=m_contained[i].size();
		for(unsigned int j=0;j<size;j++){
			rep+=mesh.positions()[m_contained[i][j]];
		}
		rep=rep/size;
		m_representatives.push_back(rep);
	}
}
