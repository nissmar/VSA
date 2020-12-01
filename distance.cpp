#include "distance.h"

//******************Distortion error******************

double triangle_area(Vector3d v1,Vector3d v2,Vector3d v3){

  double l1 = (v2-v1).norm();
  double l2 = (v3-v2).norm();
  double l3 = (v1-v3).norm();
  double p = (l1+l2+l3)/2.;
  double S = pow(p*(p-l1)*(p-l2)*(p-l3),0.5);

  return S;
};

double orthogonal_distance(Vector3d X, Vector3d N, Vector3d M){
  
  double n = N.norm();
  return fabs((N).dot(M-X))/n;

};

Vector3d triangle_normal(Vector3d v1,Vector3d v2,Vector3d v3){

  Vector3d N = (v3-v1).cross(v2-v3);
  double n = N.norm();

  return N/n;

};

Vector3d triangle_normal(Vector3i T, MatrixXd V){
  return triangle_normal(V.row(T(0)),V.row(T(1)),V.row(T(2)));
}


Vector3d triangle_center(Vector3i T, MatrixXd V){
  return (V.row(T(0)) + V.row(T(1)) + V.row(T(2))) / 3.0;
}

double distance_L_2(Vector3i T, Vector3d X, Vector3d N, MatrixXd V){

  Vector3d v1 = V.row(T(0));
  Vector3d v2 = V.row(T(1));
  Vector3d v3 = V.row(T(2));

  double area = triangle_area(v1,v2,v3);
  double d1 = orthogonal_distance(X,N,v1);
  double d2 = orthogonal_distance(X,N,v2);
  double d3 = orthogonal_distance(X,N,v3);

  return (1/6.)*area*(d1*d1 + d2*d2 + d3*d3 + d1*d2 + d1*d3 + d2*d3);

};

double distance_L_2_1(Vector3i T, Vector3d N, MatrixXd V){

  Vector3d v1 = V.row(T(0));
  Vector3d v2 = V.row(T(1));
  Vector3d v3 = V.row(T(2));
  double area = triangle_area(v1,v2,v3);
  Vector3d n = triangle_normal(v1,v2,v3);

  return area*pow((n-N).norm(),2);

};
