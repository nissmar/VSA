#include <igl/opengl/glfw/Viewer.h>
#include <queue>
#include <map>

using namespace Eigen;
using namespace std;

void tcolor(MatrixXi &Pf);
void fcolor(MatrixXd &Cf, MatrixXi Ad);

MatrixXi face_adjacency(MatrixXi F, int n);


//******************Distortion error******************
double triangle_area(Vector3d v1,Vector3d v2,Vector3d v3);
double orthogonal_distance(pair<Vector3d , Vector3d> P, Vector3d M);
Vector3d triangle_normal(Vector3d v1,Vector3d v2,Vector3d v3);

double distance_L_2 (Vector3i T, pair<Vector3d , Vector3d> P, MatrixXd V);
double distance_L_2_1(Vector3i T, pair<Vector3d , Vector3d> P, MatrixXd V);