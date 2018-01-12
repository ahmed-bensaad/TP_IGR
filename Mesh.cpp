// --------------------------------------------------------------------------
// Copyright(C) 2009-2016
// Tamy Boubekeur
//
// Permission granted to use this code only for teaching projects and
// private practice.
//
// Do not distribute this code outside the teaching assignements.
// All rights reserved.
// --------------------------------------------------------------------------

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

void Mesh::clear () {
    m_positions.clear ();
    m_normals.clear ();
    m_triangles.clear ();
    m_adj.clear();
}

void Mesh::loadOFF (const std::string & filename) {
    clear ();
	ifstream in (filename.c_str ());
    if (!in)
        exit (1);
	string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    m_positions.resize (sizeV);
    m_triangles.resize (sizeT);
    for (unsigned int i = 0; i < sizeV; i++)
        in >> m_positions[i];
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> m_triangles[i][j];
    }
    in.close ();
    centerAndScaleToUnit ();
    recomputeNormals ();
}

void Mesh::recomputeNormals () {
    m_normals.clear ();
    m_normals.resize (m_positions.size (), Vec3f (0.f, 0.f, 0.f));
    for (unsigned int i = 0; i < m_triangles.size (); i++) {
        Vec3f e01 = m_positions[m_triangles[i][1]] -  m_positions[m_triangles[i][0]];
        Vec3f e02 = m_positions[m_triangles[i][2]] -  m_positions[m_triangles[i][0]];
        Vec3f n = cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            m_normals[m_triangles[i][j]] += n;
    }
    for (unsigned int i = 0; i < m_normals.size (); i++)
        m_normals[i].normalize ();
}

void Mesh::centerAndScaleToUnit () {
    Vec3f c;
    for  (unsigned int i = 0; i < m_positions.size (); i++)
        c += m_positions[i];
    c /= m_positions.size ();
    float maxD = dist (m_positions[0], c);
    for (unsigned int i = 0; i < m_positions.size (); i++){
        float m = dist (m_positions[i], c);
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < m_positions.size (); i++)
        m_positions[i] = (m_positions[i] - c) / maxD;
}
vector<int> Mesh::getNeighbours(Vec3f vertex){
        std::vector<int> neighbours;
        for (unsigned int j=0;j<m_triangles.size();j++){
            if (m_positions[m_triangles[j][0]]==vertex){
                neighbours.push_back(m_triangles[j][1]);
                neighbours.push_back(m_triangles[j][2]);
            }
            else if (m_positions[m_triangles[j][1]]==vertex){
                neighbours.push_back(m_triangles[j][0]);
                neighbours.push_back(m_triangles[j][2]);
            }
            else if (m_positions[m_triangles[j][2]]==vertex){
                neighbours.push_back(m_triangles[j][1]);
                neighbours.push_back(m_triangles[j][0]);
            }
    }
    sort( neighbours.begin(), neighbours.end());
    neighbours.erase( unique( neighbours.begin(), neighbours.end() ), neighbours.end() );

    return neighbours;
}

void Mesh::computeAdj(){
    m_adj.resize(m_positions.size());
    for (unsigned int i=0;i<m_triangles.size();i++){
        for(unsigned int j=0;j<3;j++){
            m_adj[m_triangles[i][j]].push_back(m_triangles[i][(j+1)%3]);
            m_adj[m_triangles[i][j]].push_back(m_triangles[i][(j+2)%3]);
        }
    }
}


void Mesh::compute_cotangent_weights(){
  cotangent_weight.clear();
  cotangent_weight.resize( positions().size() );

  for(unsigned int i=0;i<m_triangles.size();i++){
    int i0 =m_triangles[i][0];
    int i1 =m_triangles[i][1];
    int i2 =m_triangles[i][2];

    Vec3f v1=m_positions[i0];
    Vec3f v2=m_positions[i1];
    Vec3f v3=m_positions[i2];

    Vec3f edge1= v2-v1;
    Vec3f edge2= v3-v2;
    Vec3f edge3= v1-v3;

    float len1=edge1.length();
    float len2=edge2.length();
    float len3=edge3.length();

    float cos1=-dot(edge1,edge3);
    float cos2=-dot(edge1,edge2);
    float cos3=-dot(edge2,edge3);

    float sin1=cross(edge1,edge3).length();
    float sin2=cross(edge1,edge2).length();
    float sin3=cross(edge2,edge3).length();


    cotangent_weight[i0][i1]+=0.5*cos3/sin3;
    cotangent_weight[i1][i0]+=0.5*cos3/sin3;

    cotangent_weight[i0][i2]+=0.5*cos2/sin2;
    cotangent_weight[i2][i0]+=0.5*cos2/sin2;

    cotangent_weight[i2][i1]+=0.5*cos1/sin1;
    cotangent_weight[i1][i2]+=0.5*cos1/sin1;


  }
}

void Mesh::compute_areas(){
  vertex_area.clear();
  vertex_area.resize( positions().size() );

  for(unsigned int i=0;i<m_triangles.size();i++){
    int i0 =m_triangles[i][0];
    int i1 =m_triangles[i][1];
    int i2 =m_triangles[i][2];

    Vec3f v1=m_positions[i0];
    Vec3f v2=m_positions[i1];
    Vec3f v3=m_positions[i2];

    Vec3f edge1= v2-v1;
    Vec3f edge2= v3-v2;

    float area=fabs(0.5*(cross(edge1,edge2)).length());
    vertex_area[i0]+=area;
    vertex_area[i1]+=area;
    vertex_area[i2]+=area;
  }

}

void Mesh::laplacianFilter(){

    for(unsigned int i=0;i<m_adj.size();i++){
        Vec3f bary=Vec3f(0.0f,0.0f,0.0f);
        for (unsigned int j=0;j<m_adj[i].size();j++){
            int index=m_adj[i][j];
            bary+=m_positions[index];
        }
        bary*=(1.f/(m_adj[i].size()));
        m_positions[i]=bary;

    }

    recomputeNormals();
}





void Mesh::GeomFilter(){

    /*Vec3f result=Vec3f(0.f,0.f,0.f);
    for(unsigned int j=0;j<m_positions.size();j++){
        Vec3f vertex=m_positions[j];
        for(unsigned int i=0;i<m_triadj[j].size();i++){
            result+=computeTriangleContrib(m_triangles[m_triadj[j][i]],vertex);
        }
        m_positions[j]+=result;
    }
    recomputeNormals();*/

    for(unsigned int i=0;i<m_positions.size();i++){
    std::map<unsigned int,float> voisins_i=cotangent_weight[i];
    for(std::map<unsigned int,float>::iterator it=voisins_i.begin();it!=voisins_i.end();++it){
      int j =it->first;
      float wij=it->second;
      m_positions[i]+=vertex_area[i]*wij*(m_positions[j]-m_positions[i])/6;
    }
  }
  recomputeNormals();

}


void Mesh::simplify(unsigned int resolution){
    float xmin= 999.f;
    float xmax= -999.f;
    float ymin= 999.f;
    float ymax= -999.f;
    float zmin= 999.f;
    float zmax= -999.f;
    std::vector<std::vector<int>> contained;
    std::vector<Vec3f> Rep;

    // Generate the big Cube and the grid
    for (unsigned int i=0;i<m_positions.size();i++){
        float x=m_positions[i][0];
        float y=m_positions[i][1];
        float z=m_positions[i][2];
        if (x<xmin){
            xmin=x;
        }
        else if (x>xmax){
            xmax=x;
        }
        if (y<ymin){
            ymin=y;
        }
        else if (y>ymax){
            ymax=y;
        }
        if (z<zmin){
            zmin=z;
        }
        else if (z>zmax){
            zmax=z;
        }
    }
    Vec3f p1=Vec3f(xmin+0.1,ymax+0.1,zmin+0.1);
    Vec3f p2=Vec3f(xmax+0.1,ymin+0.1,zmax+0.1);
    Cube container= Cube(p1,p2);
    Grid grid= Grid(container,resolution);
    grid.computeRep(*this);
    Rep = grid.getRep();
    contained = grid.getContained();

    //simplification
    for(unsigned int i=0;i<m_triangles.size();i++){
        Vec3f s1= m_positions[m_triangles[i][0]];
        Vec3f s2= m_positions[m_triangles[i][1]];
        Vec3f s3= m_positions[m_triangles[i][2]];
        unsigned int index_s1;
        unsigned int index_s2;
        unsigned int index_s3;
        for(unsigned int j=0;j<contained.size();j++){
            for (unsigned int k=0;k<contained[j].size();k++){
                if (m_positions[contained[j][k]]==s1){
                    index_s1=j;
                }
                else if(m_positions[contained[j][k]]==s2){
                    index_s2=j;
                }
                else if (m_positions[contained[j][k]]==s3){
                    index_s3=j;
                }
            }
        }
        if(index_s1==index_s2 || index_s2==index_s3 || index_s1==index_s3){
            m_triangles.erase(m_triangles.begin()+i-1);
        }
        else{
            m_triangles[i]=Triangle(index_s1,index_s2,index_s3);
        }
    }
    m_positions.clear();
    m_positions.resize(Rep.size());
    m_positions=Rep;
    recomputeNormals();
}
