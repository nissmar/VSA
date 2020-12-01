#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

//******************Distortion error******************
double triangle_area(Vector3d v1,Vector3d v2,Vector3d v3);
double orthogonal_distance(Vector3d X, Vector3d N, Vector3d M);
Vector3d triangle_normal(Vector3d v1,Vector3d v2,Vector3d v3);
Vector3d triangle_normal(Vector3i T, MatrixXd V);
Vector3d triangle_center(Vector3i T, MatrixXd V);

double distance_L_2(Vector3i T, Vector3d X, Vector3d N, MatrixXd V);
double distance_L_2_1(Vector3i T, Vector3d N, MatrixXd V);
double distance(Vector3i T, Vector3d X, Vector3d N, MatrixXd V, int norme);
