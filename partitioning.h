#include <igl/opengl/glfw/Viewer.h>
#include <queue>
using namespace Eigen;
using namespace std;

void tcolor(MatrixXi &Pf);
void fcolor(MatrixXd &Cf, MatrixXi Ad);

MatrixXi face_adjacency(MatrixXi F);
