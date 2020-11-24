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


//******************Distortion error******************
double triangle_area(Vector3d v1,Vector3d v2,Vector3d v3);
double orthogonal_distance(Vector3d X, Vector3d N, Vector3d M);
Vector3d triangle_normal(Vector3d v1,Vector3d v2,Vector3d v3);

double distance_L_2 (Vector3i T, Vector3d X, Vector3d N, MatrixXd V);
double distance_L_2_1(Vector3i T, Vector3d N, MatrixXd V);
