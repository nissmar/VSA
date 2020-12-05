#include <igl/opengl/glfw/Viewer.h>
#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

#include "distance.h"


using namespace Eigen;
using namespace std;

vector<vector<int>> anchor_points(HalfedgeDS he, MatrixXi R, MatrixXd V, MatrixXd Proxies);