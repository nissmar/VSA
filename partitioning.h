#include <igl/opengl/glfw/Viewer.h>
#include <queue>
#include <map>

using namespace Eigen;
using namespace std;

MatrixXi face_adjacency(MatrixXi F, int n);
MatrixXi uniform_proxies(int k, int n);
void tcolor(MatrixXi &R);
void fcolor(MatrixXd &Cf, MatrixXi Ad);
void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V);
void initial_partition(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad);
VectorXi find_best_triangles(MatrixXd Proxies, MatrixXd V, MatrixXi F);
void proxy_color(MatrixXi &R, MatrixXd Proxies, MatrixXd V, MatrixXi F, MatrixXi Ad);




