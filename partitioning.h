#include <igl/opengl/glfw/Viewer.h>
#include <queue>
#include <map>

using namespace Eigen;
using namespace std;

void tcolor(MatrixXi &Pf);
void fcolor(MatrixXd &Cf, MatrixXi Ad);
void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V);
void proxy_color(MatrixXi &Pf, MatrixXi Proxies, MatrixXd V, MatrixXi F, MatrixXi Ad);

MatrixXi uniform_proxies(int k, int n);
MatrixXi face_adjacency(MatrixXi F, int n);

//******************Finding new proxies******************
Vector3d g(Vector3d v1,Vector3d v2,Vector3d v3); //same notation as in appendix B
MatrixXd M(Vector3d v1,Vector3d v2,Vector3d v3); //same notation as in appendix B

Vector3d new_Xi_L_2 (MatrixXi R, int i, MatrixXi F, MatrixXd V);
Vector3d new_Ni_L_2 (MatrixXi R, int i, MatrixXi F, MatrixXd V);
Vector3d new_Xi_L_2_1 (MatrixXi R, int i, MatrixXi F, MatrixXd V);
Vector3d new_Ni_L_2_1 (MatrixXi R, int i, MatrixXi F, MatrixXd V);

MatrixXd new_proxies_L_2(MatrixXi R, MatrixXi F, MatrixXd V, int k); //k is the number of regions/proxies of the partition R
MatrixXd new_proxies_L_2_1(MatrixXi R, MatrixXi F, MatrixXd V, int k);

