#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/jet.h>
#include <iostream>
#include <ostream>

#include "partitioning.h"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F; // incidence relations between faces and edges (f columns)
MatrixXi Partition_faces; // matrix indicating the partition of each vertex
MatrixXi Ad; // face adjacency

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  cout << "pressed Key: " << key << " " << (unsigned int)key << endl;
  return false;
}


// ------------ main program ----------------
int main(int argc, char *argv[])
{
  igl::readOFF("../data/bunny.off", V, F); // Load an input mesh in OFF format

  // Face adjacency
  cout << "Computing face adjacency..." << endl;
  Ad = face_adjacency(F);
  cout << "   ...done" << endl;
  //coloring 
  // Partition_faces.setZero(F.rows(),1);
  // MatrixXd C;
  // tcolor(Partition_faces);
  // igl::jet(Partition_faces,true,C);

  // coloring 
  MatrixXd Cf;
  fcolor(Cf,Ad);
  MatrixXd C;
  igl::jet(Cf,true,C);
  
  //  print the number of mesh elements
  cout << "Vertices: " << V.rows() << endl;
  cout << "Faces:    " << F.rows() << endl;

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
  viewer.data().set_colors(C);
  viewer.launch(); // run the editor
}
