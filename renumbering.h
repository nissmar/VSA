#include <igl/opengl/glfw/Viewer.h>
#include "anchors.h"

using namespace Eigen;
using namespace std;

map<int,int> renumber(MatrixXi& F);
MatrixXd new_V(HalfedgeDS he, MatrixXd V, MatrixXd Proxies, MatrixXi R, map<int,int> index);
Vector3d projection(Vector3d z, MatrixXd Proxies, int k);