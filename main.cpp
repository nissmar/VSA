#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/jet.h>
#include <iostream>
#include <ostream>

#include "partitioning.h"
#include "distance.h"
#include "proxies.h"
#include "remeshing.h"


using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F; // incidence relations between faces and edges (f columns)
MatrixXi R; // matrix indicating the partition of each vertex
MatrixXd C; // the coloring
MatrixXd Proxies;
MatrixXi Ad; // face adjacency
int p; // number of proxies
int norme; //0 = norm L_2, 1 = norm L_2_1

void debug_regions_vides(MatrixXi R, int p){
  cout<<"Regions vides"<<endl;
  bool trouve_j;
  for (int j=0 ; j<p ; j++){
    trouve_j = false;
    for (int i=0 ; i<R.rows() ; i++){
      if (R(i,0)==j){
        trouve_j = true;
      }
    }
    if (trouve_j == false){
      cout<<j<<endl;
    }
  }
  cout<<"fin"<<endl;
};

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

void draw_anchors(igl::opengl::glfw::Viewer &viewer) {
  VectorXi anchors = anchor_points(R,Ad,F,V.rows());

  viewer.append_mesh();
  for (int i =0; i<anchors.rows(); i++) {
    viewer.data(0).add_points(V.row(anchors(i)), Eigen::RowVector3d(0, 0, 1));
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
   
    proxy_color(R, Proxies, V,  F, Ad, norme);
    Proxies = new_proxies(R, F, V, p, norme);
    
    igl::jet(R,true,C);
    viewer.data().set_colors(C);

    viewer.data(0).clear();
  }
  if (key=='4') {
    draw_anchors(viewer);
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

  //test L2
  /**Vector3d v1(0,0,1);
  Vector3d v2(1,0,1);
  Vector3d v3(0,1,1);
  MatrixXd V_test(3,3);
  V_test.row(0) = v1;
  V_test.row(1) = v2;
  V_test.row(2) = v3;
  MatrixXi F_test(1,3);
  F_test.row(0) = Vector3i(0,1,2);
  Vector3d X(0,0,0);
  Vector3d N(0,0,1);
  cout<<"distance : "<<distance_L_2_1(F_test.row(0),N,V_test)<<endl;*/

  // coloring proxies
  p = 10;
  norme = 1;
  initial_partition(p, R, V, F, Ad, norme);
  Proxies = new_proxies(R, F, V, p, norme);
  igl::jet(R,true,C);

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer

  //showing normals
  // viewer.append_mesh();
  // for (int j=0;j<F.rows();j++) {
  //   Vector3d center = triangle_center(F.row(j),V);
  //   Vector3d norm = triangle_normal(F.row(j),V);
  //   viewer.data(0).add_edges(
  //       center.transpose(),
  //       center.transpose()+norm.transpose()/10.0,
  //       Eigen::RowVector3d(1, 0, 0));
  // }

  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  draw_tangent(viewer);
  viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
  viewer.data().set_colors(C);
  viewer.launch(); // run the editor

  
}

