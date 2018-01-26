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

#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include <set>

#include "Vec3.h"
#include "Camera.h"
#include "GLProgram.h"
#include "Exception.h"
#include "LightSource.h"
#include "Ray.h"
#include "Vec4.h"

using namespace std;

static const unsigned int DEFAULT_SCREENWIDTH = 1024;
static const unsigned int DEFAULT_SCREENHEIGHT = 768;
static const string DEFAULT_MESH_FILE ("models/sphere.off");

static const string appTitle ("Informatique Graphique & Realite Virtuelle - Travaux Pratiques - Algorithmes de Rendu");
static const string myName ("Ahmed Ben Saad");
static GLint window;
static unsigned int FPS = 0;
static bool fullScreen = false;
static Camera camera;
static Mesh mesh;
GLProgram * glProgram;
GLProgram * toonProgram;
GLProgram * toon_outlines_Program;
GLProgram * bvhProgram;

static std::vector<Vec4f> colorResponses;		// Cached per-vertex color response, updated at each frame
static std::vector<LightSource> lightSources;

static BVH * big_bvh;
static unsigned int currentDeep = 0;

std::vector<Vec3f> BVH::bvh_positions;
std::vector<unsigned int> BVH::bvh_indices;
unsigned int BVH::deep_count = 0;


static int nb_Samples = 20;		 //number of sample of rays in one point when calculating AO
static float radius_AO = 0.6f;   //radius of rays when calculating AO

static bool cartoon_mode = false;   //cartoon mode
static bool perVertexShadow = false;  //Shadow mode (true by default)
static bool renderShadowOnlyInInit = true;  //Render shadow only in the beginning of the program or in every frame ( calls of function renderScene() )
static bool perVertexAO = false;     //Ambient occlusion mode (true by default)
static bool ggx = true;					//Cook-Torrance micro facet BRDF / GGX micro facet BRDF
static bool filter=false;
void printUsage () {
	std::cerr << std::endl<< appTitle << std::endl
         << "Author: " << myName << std::endl << std::endl
         << "Usage: ./main [<file.off>]" << std::endl
         << "Commands:" << std::endl
         << "------------------" << std::endl
         << " ?: Print help" << std::endl
		 << " w: Toggle wireframe mode" << std::endl
         << " <drag>+<left button>: rotate model" << std::endl
         << " <drag>+<right button>: move model" << std::endl
         << " <drag>+<middle button>: zoom" << std::endl
				 << " <f>: full screen mode"<< std::endl
				 << " <t>: active cartoon mode"<< std::endl
        << " <b>: display bounding boxes"<< std::endl
					<< " <c>: Geometric filter"<< std::endl
					<< " <x>: Laplacian filter"<< std::endl
        << " <v>: Grid simplification 16x16x16"<< std::endl
        << " <n>: Grid simplification 32x32x32"<< std::endl
        << " <o>: Grid simplification 64x64x64"<< std::endl

         << " <q>, <esc>: Quit" << std::endl << std::endl;
}
void polar2Cartesian (float phi, float theta, float r, float & x, float & y, float & z) {
    x = r * sin (theta) * cos (phi);
    y = r * sin (theta) * sin (phi);
    z = r * cos (theta);
  }

void plan(){
    std::vector<Vec3f> & positions = mesh.positions();
    std::vector<Triangle> & triangles = mesh.triangles();
    unsigned int N = positions.size();
    unsigned int n = 40;
    for (unsigned int i = 0 ; i < n ; i++){
    for (unsigned int j = 0 ; j < n ; j++){
        double ui = (double)(i) / (double)(n-1);
        double uj = (double)(j) / (double)(n-1);
    positions.push_back(Vec3f(-2. + 4*uj,-1.0,-2. + 4*ui));
         }
    }
    for (unsigned int j = 0; j < n-1 ; j++){
        for (unsigned int i = 0; i < n-1; i++){
            triangles.push_back(Triangle(N+j+n*i,N+j+n*(i+1),N+j+1+n*i)) ;
            triangles.push_back(Triangle(N+j+n*(i+1),N+j+1+n*(i+1),N+j+1+n*i)) ;
        }
}
mesh.recomputeNormals();
}
void computePerVertexAO (unsigned int numOfSamples, float radius, const Mesh& mesh){
    srand((unsigned)time(0));
    unsigned int mesh_size = mesh.positions().size();
    for(unsigned int i = 0; i < mesh_size; i++){
        Vec3f p = mesh.positions()[i];
        Vec3f n = mesh.normals()[i];
        float mark = 0.0f;
        unsigned int k = 0;
        while(true){
            float random_variable1 = (float)rand() / (RAND_MAX);
            float random_variable2 = (float)rand() / (RAND_MAX);
            float random_variable3 = (float)rand() / (RAND_MAX);
            int sign1 = rand() % 2;
            int sign2 = rand() % 2;
            int sign3 = rand() % 2;
            random_variable1 = (sign1 == 0) ? random_variable1 : - random_variable1;
            random_variable2 = (sign2 == 0) ? random_variable2 : - random_variable2;
            random_variable3 = (sign3 == 0) ? random_variable3 : - random_variable3;
            Vec3f direction = Vec3f(random_variable1, random_variable2, random_variable3);
            if(  length(direction) > 1.0 ) continue;
            if(  ( dot(direction, n) / ( length(direction) * length(n) ) ) < 0.7f ) continue;
            k++;
            direction.normalize();
            direction = direction * radius;
            Ray ray(p, direction);
            if( ray.isIntersected(mesh, big_bvh) ) mark += 1.0f;
            if(k == numOfSamples) break;
        }
        colorResponses[i][0] = (float)(mark / numOfSamples);
    }
}

void computePerVertexShadow(const Mesh& mesh){
    unsigned int ls_size = 1;
    unsigned int mesh_size = mesh.positions().size();
    for(unsigned int k = 0; k < ls_size; k++){
        Vec3f lightPosition = lightSources[k].getPosition();
        for(unsigned int i = 0; i < mesh_size; i++){
            Vec3f origin = mesh.positions()[i];
        Vec3f direction = lightPosition - origin;
            Ray ray(origin + direction * 0.00001f, direction);
            if( ray.isIntersected(mesh, big_bvh) ){
                colorResponses[i][3] = 0.0f;
            }else{
                colorResponses[i][3] = 1.0f;
            }
        }
    }
}
void init (const char * modelFilename) {
    glewExperimental = GL_TRUE;
    glewInit (); // init glew, which takes in charges the modern OpenGL calls (v>1.2, shaders, etc)
    glCullFace (GL_BACK);     // Specifies the faces to cull (here the ones pointing away from the camera)
    glEnable (GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
    glDepthFunc (GL_LESS); // Specify the depth test for the z-buffer
    glEnable (GL_DEPTH_TEST); // Enable the z-buffer in the rasterization
    glEnableClientState (GL_VERTEX_ARRAY);
    glEnableClientState (GL_NORMAL_ARRAY);
    glEnableClientState (GL_COLOR_ARRAY);
    glEnable (GL_NORMALIZE);
	glLineWidth (2.0); // Set the width of edges in GL_LINE polygon mode
    glClearColor (0.0f, 0.0f, 0.0f, 1.0f); // Background color

	mesh.loadOFF (modelFilename);
	plan();			//adds a plane
    mesh.computeAdj();
	mesh.compute_areas();
	mesh.compute_cotangent_weights();
    colorResponses.resize (mesh.positions().size());
    camera.resize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);

    try {
        glProgram = GLProgram::genVFProgram ("Simple GL Program", "shader.vert", "shader.frag"); // Load and compile pair of shaders
				toonProgram = GLProgram::genVFProgram ("Cartoon GL Program", "shader_cartoon.vert", "shader_cartoon.frag");
				toon_outlines_Program = GLProgram::genVFProgram ("Outline GL Program", "shader_outline_cartoon.vert", "shader_outline_cartoon.frag");
				bvhProgram = GLProgram::genVFProgram ("BVH GL Program", "shader_bvh.vert", "shader_bvh.frag");
				glProgram->use (); // Activate the shader program
    } catch (Exception & e) {
        cerr << e.msg () << endl;
    }

        unsigned int deep_count1 = 0;
		big_bvh = BVH::buildBVH( mesh.triangles(), mesh, deep_count1);

		lightSources.resize(8);
		lightSources[0] = LightSource(Vec3f(0.0f, 3.0f, -3.0f), Vec3f(1.0f, 0.5f, 0.2f));
		lightSources[0].activeLightSource();
        lightSources[1] = LightSource(Vec3f(-3.0f, 4.0f, 0.0f), Vec3f(0.0f, 1.0f, 1.0f));
        lightSources[1].activeLightSource();
		if(renderShadowOnlyInInit == true){
			if(perVertexAO == true) computePerVertexAO(nb_Samples, radius_AO,mesh);
			if(perVertexShadow == true) computePerVertexShadow(mesh);
		}

}

/*
void updatePerVertexColorResponse () {

    computePerVertexShadow ();
    computePerVertexAO(20,0.05f);

    for (unsigned int i = 0; i < colorResponses.size (); i++){
        Vec3f light0_dir=light0.Position-mesh.positions()[i];
        float d0= light0_dir.squaredLength();

        Vec3f light1_dir=light1.Position-mesh.positions()[i];
        float d1= light1_dir.squaredLength();


        Vec3f light2_dir=light2.Position-mesh.positions()[i];
        float d2= light2_dir.squaredLength();


        Vec3f cam;
        camera.getPos(cam);
        Vec3f emit_dir=(cam-mesh.positions()[i]);
        emit_dir.normalize();
        mesh.normals()[i].normalize();


        Vec3f wh0=(light0_dir+emit_dir);
        wh0.normalize();
        Vec3f wh1=(light1_dir+emit_dir);
        wh1.normalize();
        Vec3f wh2=(light2_dir+emit_dir);
        wh2.normalize();


        //colorResponses[i]=((float)(max<float>(dot(mesh.normals()[i],light0_dir),0.f)*(1.f/M_PI))*light0.Color)*(light0.intensity/d0) +((float)(max<float>(dot(mesh.normals()[i],light1_dir),0.f)*(1.f/M_PI))*light1.Color)*(light1.intensity/d1)+((float)(max<float>(dot(mesh.normals()[i],light2_dir),0.f)*(1.f/M_PI))*light2.Color)*(light2.intensity/d2) ;



        Vec3f colorResponse0=(max<float>(dot(mesh.normals()[i],light0_dir),0.f)*light0.Color)*(light0.intensity/d0)*pow((dot(mesh.normals()[i],wh0)),3);
        Vec3f colorResponse1=(max<float>(dot(mesh.normals()[i],light1_dir),0.f)*light1.Color)*(light1.intensity/d1)*pow((dot(mesh.normals()[i],wh1)),3);
        Vec3f colorResponse2=(max<float>(dot(mesh.normals()[i],light2_dir),0.f)*light2.Color)*(light2.intensity/d2)*pow((dot(mesh.normals()[i],wh2)),3);
        colorResponses[i]=colorResponse2+colorResponse1+colorResponse0;

        //w0 vers la cam = emit_dir.n
        //wh=wi+w0/norme
        //wi vers la lumiere=light_dir.n

        float n_wh1=max(0.f,dot(wh1,mesh.normals()[i]));
        float n_wh2=max(0.f,dot(wh2,mesh.normals()[i]));
        float n_wi1=max(0.f,dot(mesh.normals()[i],light1_dir));
        float n_wi2=max(0.f,dot(mesh.normals()[i],light2_dir));








       float D1=(1.f/(M_PI*pow(alpha,2)*pow(n_wh1,4)))*exp((pow(n_wh1,2)-1)/(pow(alpha,2)*pow(n_wh1,2)));
        float F1=Fi+(1-Fi)*pow(1-max(0.f,dot(emit_dir,light1_dir)),5);
        float G1=min(1.f,min((2*n_wi1*n_wh1)/dot(emit_dir,wh1),(2*n_w*n_wh1)/dot(emit_dir,wh1)));
        float D1_GGX=pow(alpha,2)/(M_PI*pow((1+(pow(alpha,2)-1)*pow(n_wh1,2)),2));
        float G1_GGX=(4*n_w*n_wi1)/((n_w+sqrt(pow(alpha,2)+(1-pow(alpha,2))*pow(n_w,2)))*(n_wi1+sqrt(pow(alpha,2)+(1-pow(alpha,2))*pow(n_wi1,2))));
        float MF1=(D1_GGX*F1*G1_GGX)/(4*n_wi1*n_w);


        float D2=(1.f/(M_PI*pow(alpha,2)*pow(n_wh2,4)))*exp((pow(n_wh2,2)-1)/(pow(alpha,2)*pow(n_wh2,2)));
        float F2=Fi+(1-Fi)*pow(1-max(0.f,dot(emit_dir,light2_dir)),5);
        float G2=min(1.f,min((2*n_wi2*n_wh2)/dot(emit_dir,wh2),(2*n_w*n_wh2)/dot(emit_dir,wh2)));
        float D2_GGX=pow(alpha,2)/(M_PI*pow((1+(pow(alpha,2)-1)*pow(n_wh2,2)),2));
        float G2_GGX=(4*n_w*n_wi2)/((n_w+sqrt(pow(alpha,2)+(1-pow(alpha,2))*pow(n_w,2)))*(n_wi2+sqrt(pow(alpha,2)+(1-pow(alpha,2))*pow(n_wi2,2))));
        float MF2=(D2_GGX*F2*G2_GGX)/(4*n_wi2*n_w);

        //+(max<float>(dot(mesh.normals()[i],light1_dir),0.f)*light1.Color)*MF1*(light1.intensity/d1)+(max<float>(dot(mesh.normals()[i],light2_dir),0.f)*light2.Color)*MF2*(light2.intensity/d2);


*/

void renderScene () {
    //updatePerVertexColorResponse ();

		GLfloat lightPos[3 * 8];
		GLfloat lightCol[3 * 8];
		GLint nb_light_active = 0;
		for(unsigned int i = 0; i < 8; i++){
			LightSource lightSource = lightSources[i];
			if(lightSource.isActive() == true){
				lightPos[3 * i + 0] = lightSource.getPosition()[0];
				lightPos[3 * i + 1] = lightSource.getPosition()[1];
				lightPos[3 * i + 2] = lightSource.getPosition()[2];
				lightCol[3 * i + 0] = lightSource.getColor()[0];
				lightCol[3 * i + 1] = lightSource.getColor()[1];
				lightCol[3 * i + 2] = lightSource.getColor()[2];
				nb_light_active++;
			}
		}

		if( cartoon_mode == false ){
			glProgram->use ();
			GLint variableLocationPos = glProgram->getUniformLocation("lightPositions");
			GLint variableLocationCol = glProgram->getUniformLocation("lightColors");
			glUniform3fv (variableLocationPos, nb_light_active, lightPos);
			glUniform3fv (variableLocationCol, nb_light_active, lightCol);
			glProgram->setUniform1i("numberOfLightActive", nb_light_active);
			glProgram->setUniform1i("ggx", ggx);
			glProgram->setUniform1i("perVertexShadow", perVertexShadow);
			glProgram->setUniform1i("perVertexAO", perVertexAO);
		}else{
			glCullFace (GL_FRONT);
			//toon_outlines_Program->use();
			glVertexPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(mesh.positions()[0])));
	    glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)&(mesh.normals()[0]));
			glColorPointer (4, GL_FLOAT, sizeof (Vec4f), (GLvoid*)(&(colorResponses[0])));
			glDrawElements (GL_TRIANGLES, 3*mesh.triangles().size(), GL_UNSIGNED_INT, (GLvoid*)((&mesh.triangles()[0])));

			glCullFace (GL_BACK);
			toonProgram->use();
			GLint variableLocationPos = toonProgram->getUniformLocation("lightPositions");
			glUniform3fv (variableLocationPos, nb_light_active, lightPos);
			toonProgram->setUniform1i("numberOfLightActive", nb_light_active);

			toonProgram->setUniform1i("perVertexShadow", perVertexShadow);
			toonProgram->setUniform1i("perVertexAO", perVertexAO);
		}

		if(renderShadowOnlyInInit == false){
			if(perVertexShadow == true) computePerVertexShadow(mesh);
		}

    glVertexPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(mesh.positions()[0])));
    glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)&(mesh.normals()[0]));
		glColorPointer (4, GL_FLOAT, sizeof (Vec4f), (GLvoid*)(&(colorResponses[0])));
    glDrawElements (GL_TRIANGLES, 3*mesh.triangles().size(), GL_UNSIGNED_INT, (GLvoid*)((&mesh.triangles()[0])));

		bvhProgram->use();
		glVertexPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(BVH::bvh_positions[0])));
		glDrawElements (GL_LINES, BVH::bvh_indices.size(), GL_UNSIGNED_INT, (GLvoid*)((&BVH::bvh_indices[0])));
}

void reshape(int w, int h) {
    camera.resize (w, h);
}

void display () {
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    renderScene ();
    glFlush ();
    glutSwapBuffers ();
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen) {
            glutReshapeWindow (camera.getScreenWidth (), camera.getScreenHeight ());
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;
    case 'q':
    case 27:
        exit (0);
        break;
    case 'w':
        GLint mode[2];
		glGetIntegerv (GL_POLYGON_MODE, mode);
		glPolygonMode (GL_FRONT_AND_BACK, mode[1] ==  GL_FILL ? GL_LINE : GL_FILL);
        break;
        break;
		case 't':
				cartoon_mode = ! cartoon_mode;
				break;
        case 'b':
                if( currentDeep < 15 ){
                    currentDeep++;
                    big_bvh->drawBVH(currentDeep);
                }
                else{
                    currentDeep=1;
                }
                break;
	   case 'x':{
                    mesh.laplacianFilter();
                break;
       }
       case 'c':{
        mesh.GeomFilter();
        break;
       }
	   case 'v':{
            mesh.simplify(16);
            mesh.setmesh();
            renderScene();
            break;

       }
        case 'n':{
            mesh.simplify(32);
            mesh.setmesh();
            renderScene();
            break;

       }
        case 'o':{
            mesh.simplify(64);
            mesh.setmesh();
            renderScene();
            break;

       }

    default:
        printUsage ();
        break;
    }
}



void mouse (int button, int state, int x, int y) {
    camera.handleMouseClickEvent (button, state, x, y);
}

void motion (int x, int y) {
    camera.handleMouseMoveEvent (x, y);
}

void idle () {
    static float lastTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char winTitle [128];
        unsigned int numOfTriangles = mesh.triangles ().size ();
        sprintf (winTitle, "Number Of Triangles: %d - FPS: %d", numOfTriangles, FPS);
        string title = appTitle + " - By " + myName  + " - " + winTitle;
        glutSetWindowTitle (title.c_str ());
        lastTime = currentTime;
    }
    glutPostRedisplay ();
}



int main (int argc, char ** argv) {
    if (argc > 2) {
        printUsage ();
        exit (1);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    window = glutCreateWindow (appTitle.c_str ());
    init (argc == 2 ? argv[1] : DEFAULT_MESH_FILE.c_str ());
    glutIdleFunc (idle);
    glutReshapeFunc (reshape);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    printUsage ();
    glutMainLoop ();
    return 0;
}
