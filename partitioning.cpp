#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;


void color(MatrixXi &Pf) {
  int n = Pf.rows();
  for (int i=0; i<n;i++) {
    Pf(i,0) = i;
  }
}
