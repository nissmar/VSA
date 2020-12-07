#include <igl/opengl/glfw/Viewer.h>
#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

#include <queue>
#include "partitioning.h"


using namespace Eigen;
using namespace std;

MatrixXi color_region (MatrixXi R, int region, vector<vector<int>> anchors, MatrixXd V, HalfedgeDS he);
pair<MatrixXi,MatrixXi> triangulation (MatrixXi& R, vector<vector<int>> anchors, MatrixXd V, MatrixXi F, HalfedgeDS he);