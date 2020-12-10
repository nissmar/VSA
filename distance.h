#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

//******************Distortion error******************
double orthogonal_distance(Vector3d X, Vector3d N, Vector3d M);
Vector3d get_center(int i);
Vector3d get_normal(int i);
double get_area(int i);


double distance(int i, Vector3d X, Vector3d N, MatrixXd V, int norme);

double global_distortion_error(MatrixXi R, MatrixXd Proxies, MatrixXd V, MatrixXi F, int norme);

double distance_projection(MatrixXd V, MatrixXd Proxies, int anchor1, int anchor2, int v, int r1, int r2);

void initialize_normals_areas(MatrixXi F, MatrixXd V);