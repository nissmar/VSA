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
    draw_anchors(viewer);

    cout<<"\nregion 23"<<endl;
    for (int f=0 ; f<F.rows() ; f++){
      if (R(f,0)==23){
        cout<<"face "<<f<<"("<<F(f,0)<<","<<F(f,1)<<","<<F(f,2)<<")"<<endl;
      }
    }

    vector<vector<int>> anchors = anchor_points(*he, R, V, Proxies);
    //for (int i=0 ; i<anchors.size() ; i++){
      cout<<"\nregion "<<23<<" "<<anchors[23].size()<<" anchors :"<<endl;
      for (int j=0 ; j<anchors[23].size() ; j++){
        cout<<anchors[23][j]<<endl;
      }
    //}
    //triangulate_region(R,98,anchors,V,F,*he);
    //pair<MatrixXi,MatrixXi> new_F_and_R = triangulation(R,anchors,V,F,*he);
    //F = new_F_and_R.first;
    //cout<<"\nF\n"<<F<<endl;
    //R = new_F_and_R.second;
    //cout<<"\nR\n"<<R<<endl;

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


  //petit carré plan simple (5x5) pour tester la triangulation
  int n = 7;

  MatrixXd V_bis(n*n,3);
  for (int i=0 ; i<n ; i++){
    for (int j=0 ; j<n ; j++){
      V_bis.row(j+n*i) = Vector3d(i,j,0);
    }
  }

  MatrixXi F_bis((n-1)*(n-1)*2,3);
  for (int i=0 ; i<n-1 ; i++){
    for (int j=0 ; j<n-1 ; j++){
      F_bis.row(j+(n-1)*i) = Vector3i(j+n*i,j+n*(i+1),j+1+n*i);
      F_bis.row(j+(n-1)*i+(n-1)*(n-1)) = Vector3i(j+n*(i+1),j+1+n*(i+1),j+1+n*i);
    }
  }

  HalfedgeDS* he_bis;
  HalfedgeBuilder* builder_bis = new HalfedgeBuilder();  
  HalfedgeDS he2_bis = builder_bis->createMesh(V_bis.rows(), F_bis); 
  he_bis = &he2_bis;

  MatrixXi R_bis = MatrixXi::Ones((n-1)*(n-1)*2,1);
  for (int j=0 ; j<n-1 ; j++){
    R_bis(j,0) = 0;
    R_bis(j+(n-1)*(n-1),0) = 0;
    R_bis(j+(n-1)*(n-2),0) = 0;
    R_bis(j+(n-1)*(n-2)+(n-1)*(n-1),0) = 0;
  }
  for (int i=0 ; i<n-1 ; i++){
    R_bis((n-1)*i,0) = 0;
    R_bis((n-1)*i+(n-1)*(n-1),0) = 0;
    R_bis(n-2+(n-1)*i,0) = 0;
    R_bis(n-2+(n-1)*i+(n-1)*(n-1),0) = 0;
  }

  vector<int> anchors_0{0,n*(n-1),n*n-1,n-1};
  vector<int> anchors_1{n+1,(n-1)*(n-1),(n-1)*(n-1)+n-3,2*n-2};
  vector<vector<int>> anchors{anchors_0,anchors_1};
  
  /**int edge = 40;
  cout<<"edge "<<edge<<endl;
  cout<< "first pointe vers " <<he_bis->getTarget(find_first(*he_bis,edge,1,R_bis)) <<endl;
  cout<< "first pointe contre " <<he_bis->getTarget(he_bis->getOpposite(find_first(*he_bis,edge,1,R_bis))) <<endl;
  cout<< "second pointe vers " <<he_bis->getTarget(find_second(*he_bis,edge,1,R_bis)) <<endl;
  cout<< "second pointe contre " <<he_bis->getTarget(he_bis->getOpposite(find_second(*he_bis,edge,1,R_bis))) <<endl;
  cout<<"suite first "<<he_bis->getTarget(find_next_first(*he_bis,find_first(*he_bis,edge,1,R_bis),1,R_bis))<<endl;
  cout<<"suite second "<<he_bis->getTarget(find_next_second(*he_bis,find_second(*he_bis,edge,1,R_bis),1,R_bis))<<endl;*/
  
  cout<<"\ntriangulation"<<endl;
  vector<Vector3i> triangles = triangulate_region(R_bis,1,anchors,V_bis,F_bis,*he_bis);
  for (int i=0 ; i<triangles.size() ; i++){
    cout<<triangles[i]<<"\n"<<endl;
  }
  cout<<"\n"<<endl;

 /** vector<int> neighb = find_interior_neighbors(*he_bis,8,1,R_bis);
  cout<<neighb.size()<<" neighbors"<<endl;
  for (int i=0 ; i<neighb.size() ; i++){
    cout<<he_bis->getTarget(neighb[i])<<endl;
  }*/
  

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

