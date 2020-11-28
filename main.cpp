#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/jet.h>
#include <iostream>
#include <ostream>

#include "partitioning.h"
#include "distance.h"
#include "proxies.h"


using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F; // incidence relations between faces and edges (f columns)
MatrixXi R; // matrix indicating the partition of each vertex
MatrixXd C; // the coloring
MatrixXd Proxies;
MatrixXi Ad; // face adjacency
int p; // number of proxies

void draw_tangent(igl::opengl::glfw::Viewer &viewer) {
   for (int i =0; i<p;i++) {
    viewer.append_mesh();
    viewer.data(0).add_points(Proxies.row(i), Eigen::RowVector3d(1, 0, 0));
    viewer.data(0).add_edges(
        Proxies.row(i),
        Proxies.row(i) + Proxies.row(i+p)/10.0,
        Eigen::RowVector3d(1, 0, 0));
  }
}
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  cout << "pressed Key: " << key << " " << (unsigned int)key << endl;
  if (key=='1') {
    viewer.data().clear();
  }
  if (key=='2') {
    viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
    viewer.data().set_colors(C);
  }
  if (key=='3') {
    proxy_color(R, Proxies, V,  F, Ad);
    Proxies = new_proxies_L_2_1(R, F, V, p);
    cout << Proxies << endl;
    igl::jet(R,true,C);
    viewer.data().set_colors(C);

    viewer.data(0).clear();
    // draw_tangent(viewer);
  }
  return false;
}



// ------------ main program ----------------
int main(int argc, char *argv[])
{
  igl::readOFF("../data/bunny.off", V, F); // Load an input mesh in OFF format

  
  //  print the number of mesh elements
  cout << "Vertices: " << V.rows() << endl;
  cout << "Faces:    " << F.rows() << endl;

  // Face adjacency
  cout << "Computing face adjacency..." << endl;
  Ad = face_adjacency(F,V.rows());
  cout << "   ...done" << endl;

  //coloring 
  // Partition_faces.setZero(F.rows(),1);
  // MatrixXd C;
  // tcolor(Partition_faces);
  // igl::jet(Partition_faces,true,C);

  // coloring adjacency 
  // MatrixXd Cf;
  // fcolor(Cf,Ad);
  // MatrixXd C;
  // igl::jet(Cf,true,C);

  // coloring distance 
  // MatrixXd Cf;
  // distance_color(Cf,F,V);
  // MatrixXd C;
  // igl::jet(Cf,true,C);

  // coloring proxies
  p = 30;
  initial_partition(p, R, V,  F, Ad);
  Proxies = new_proxies_L_2(R, F, V, p);
  igl::jet(R,true,C);

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  draw_tangent(viewer);
  viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
  viewer.data().set_colors(C);
  viewer.launch(); // run the editor
}
