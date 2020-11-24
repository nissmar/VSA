#include "partitioning.h"

MatrixXi face_adjacency(MatrixXi F, int n) { // n: number of vertices
  // face adjacency matrix, O(m log(m)) performances
 
  int m = F.rows();
  MatrixXi Ad;
  Ad.setZero(m,3);
  map<int,int> edge;
  
  for (int i=0;i<m;i++) {
    int r[3] = {F(i,0),F(i,1),F(i,2)};
    sort(r, r + 3);
    // edge mapping: e -> e.min + n*e.max
    int e[3] = {r[0] + n*r[1],r[0] + n*r[2],r[1] + n*r[2]}; // order edges by vertices
    for (int l=0; l<3; l++){
      if (edge.count(e[l])!=0) { // if edge contains(e[l])
        int j = edge[e[l]];
        Ad(i,l) = j; // faces i and j are neighbourgs

        int r2[3] = {F(j,0),F(j,1),F(j,2)};
        sort(r2, r2 + 3);
        int e2[3] = {r2[0] + n*r2[1],r2[0] + n*r2[2],r2[1] + n*r2[2]};
        int k = (e[l]==e2[1]) + 2*(e[l]==e2[2]); // check corresponding index for j
        Ad(j,k) = i; // faces i and j are neighbourgs
      }
      else {
        edge[e[l]] = i;
      }
    }
  }
  return Ad;
}

void tcolor(MatrixXi &Pf) {
  int n = Pf.rows();
  for (int i=0; i<n;i++) {
    Pf(i,0) = i;
  }
}

void fcolor(MatrixXd &Cf, MatrixXi Ad) {
  Cf.setZero(Ad.rows(), 1);
  double color = 1.0;  
  priority_queue<pair<double, int>> q;
  pair<double, int> face;
  q.push(make_pair(color,0));
  while (q.size()!=0) {
    face = q.top();
    q.pop();
    if (Cf(face.second,0)==0) {
      Cf(face.second,0) = face.first;
      color = face.first * 0.95;
      q.push(make_pair(color,Ad(face.second,0)));
      q.push(make_pair(color,Ad(face.second,1)));
      q.push(make_pair(color,Ad(face.second,2)));
    }
  }
  
}

//******************Distortion error******************

double triangle_area(Vector3d v1,Vector3d v2,Vector3d v3){

  double l1 = (v2-v1).norm();
  double l2 = (v3-v2).norm();
  double l3 = (v1-v3).norm();
  double p = (l1+l2+l3)/2.;
  double S = pow(p*(p-l1)*(p-l2)*(p-l3),0.5);

  return S;
};

double orthogonal_distance(pair<Vector3d , Vector3d> P, Vector3d M){

  Vector3d X = P.first;
  Vector3d N = P.second;

  return fabs((N).dot(M-X))/N.norm();

};

Vector3d triangle_normal(Vector3d v1,Vector3d v2,Vector3d v3){

  Vector3d N = (v3-v1).cross(v2-v3);
  double n = N.norm();

  return N/n;

};

double distance_L_2 (Vector3i T, pair<Vector3d , Vector3d> P, MatrixXd V){

  Vector3d v1 = V.row(T(0));
  Vector3d v2 = V.row(T(1));
  Vector3d v3 = V.row(T(2));

  double area = triangle_area(v1,v2,v3);
  double d1 = orthogonal_distance(P,v1);
  double d2 = orthogonal_distance(P,v2);
  double d3 = orthogonal_distance(P,v3);

  return (1/6.)*area*(d1*d1 + d2*d2 + d3*d3 + d1*d2 + d1*d3 + d2*d3);

};

double distance_L_2_1(Vector3i T, pair<Vector3d , Vector3d> P, MatrixXd V){

  Vector3d v1 = V.row(T(0));
  Vector3d v2 = V.row(T(1));
  Vector3d v3 = V.row(T(2));

  Vector3d N = P.second;

  double area = triangle_area(v1,v2,v3);
  Vector3d n = triangle_normal(v1,v2,v3);

  return area*pow((n-N).norm(),2);

};
