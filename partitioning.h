#include <igl/opengl/glfw/Viewer.h>
#include <queue>
#include <map>

using namespace Eigen;
using namespace std;

void tcolor(MatrixXi &R);
void fcolor(MatrixXd &Cf, MatrixXi Ad);
void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V);
void initial_partition(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad);
void proxy_color(MatrixXi &R, MatrixXd Proxies, MatrixXd V, MatrixXi F, MatrixXi Ad);

MatrixXi face_adjacency(MatrixXi F, int n);
