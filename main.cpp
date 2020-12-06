#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/jet.h>
#include <iostream>
#include <ostream>

#include "HalfedgeBuilder.cpp"

#include "partitioning.h"
#include "distance.h"
#include "proxies.h"
#include "anchors.h"
#include "triangulation.h"
#include "renumbering.h"

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
HalfedgeDS* he;
int iterations;
vector<pair<int,double>> global_error_points; //contains the global_distortion_error according to the number of iterations
double error;
double precedent_error;

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
  vector<vector<int>> anchors = anchor_points(*he, R, V, Proxies);
  viewer.append_mesh();
  for(size_t i = 0; i < anchors.size(); i++) {
    for(size_t j = 0; j < anchors[i].size(); j++) {
      viewer.data(0).add_points(V.row(anchors[i][j]), Eigen::RowVector3d(i%3, (i+1)%3,  (i+2)%3));
    }
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
    // debug_regions_vides(R,p);
    proxy_color(R, Proxies, V,  F, Ad, norme);
    Proxies = new_proxies(R, F, V, p, norme);
    iterations += 1;
    double error = global_distortion_error(R,Proxies,V,F,norme);
    cout<<"Global Error : "<<error<<endl;
    global_error_points.push_back(make_pair(iterations,error));

    igl::jet(R,true,C);
    viewer.data(0).set_colors(C);
  }
  if (key=='4') {
    draw_anchors(viewer);
  }
  if (key=='5') {

    vector<vector<int>> anchors = anchor_points(*he, R, V, Proxies);
    pair<MatrixXi,MatrixXi> new_F_and_R = triangulation(R,anchors,V,F,*he);
    F = new_F_and_R.first;
    R = new_F_and_R.second;

    viewer.data().clear();
    igl::jet(R,true,C);
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
  }
  if (key == 'S' || (unsigned int)key == 83){
    while (fabs(error - precedent_error)>0.0001){

    proxy_color(R, Proxies, V,  F, Ad, norme);
    Proxies = new_proxies(R, F, V, p, norme);
    iterations += 1;
    precedent_error = error;
    error = global_distortion_error(R,Proxies,V,F,norme);
    global_error_points.push_back(make_pair(iterations,error));

    igl::jet(R,true,C);
    viewer.data(0).set_colors(C);

    }
  }
  return false;
}



// ------------ main program ----------------
int main(int argc, char *argv[])
{
  igl::readOFF("../data/bunny.off", V, F); // Load an input mesh in OFF format
  HalfedgeBuilder* builder=new HalfedgeBuilder();  
  HalfedgeDS he2 = builder->createMesh(V.rows(), F); 
  he = &he2;
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
  p = 100;
  norme = 1;
  initial_partition(p, R, V, F, Ad, norme);
  Proxies = new_proxies(R, F, V, p, norme);
  iterations = 1;
  error = global_distortion_error(R,Proxies,V,F,norme);
  precedent_error = error - 1 ; 
  global_error_points.push_back(make_pair(iterations,error));
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
  viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
  viewer.data().set_colors(C);
  viewer.launch(); // run the editor


  cout<<"\n_______Erreurs par itération_______\n"<<endl;
  pair<int,double> item;
  for (int i=0 ; i<global_error_points.size() ; i++){
    item = global_error_points[i];
    cout<<"( "<<item.first<<" , "<<item.second<<" )"<<endl;
  }
}

