#include <igl/opengl/glfw/Viewer.h>
#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

#include <queue>
#include "partitioning.h"


using namespace Eigen;
using namespace std;

int find_first(HalfedgeDS he, int v, int r, MatrixXi R);
int find_second(HalfedgeDS he, int v, int r, MatrixXi R);
int find_next_first(HalfedgeDS he, int edge, int r, MatrixXi R);
int find_next_second(HalfedgeDS he, int edge, int r, MatrixXi R);
vector<int> find_interior_neighbors(HalfedgeDS he, int v, int r, MatrixXi R);

double get_length_edge(HalfedgeDS he, int edge, MatrixXd V);

int find_englobing_region(int region, MatrixXi R, HalfedgeDS he, MatrixXi F, MatrixXd V);
MatrixXi update_englobing_region(MatrixXi R_bis,int region,int englobing_region);

MatrixXi color_region (MatrixXi R, int region, vector<vector<int>> anchors, MatrixXd V, HalfedgeDS he);
vector<Vector3i> triangulate_region (MatrixXi R, int region, vector<vector<int>> anchors, MatrixXd V, MatrixXi F, HalfedgeDS he);
pair<MatrixXi,MatrixXi> triangulation (MatrixXi R, vector<vector<int>> anchors, MatrixXd V, MatrixXi F, HalfedgeDS he);