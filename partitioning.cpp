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

//******************Finding new proxies******************

Vector3d g(Vector3d v1,Vector3d v2,Vector3d v3){

  return 0.5*(v1+v2+v3);

};

MatrixXd M(Vector3d v1,Vector3d v2,Vector3d v3){

  MatrixXd M = MatrixXd::Zero(3,3);
  M.row(0) = v2-v1;
  M.row(1) = v3-v1;

  return M;

}; 

Vector3d new_Xi_L_2 (MatrixXi R, int i, MatrixXi F, MatrixXd V){

  Vector3d Xi(0.,0.,0.);
  double w = 0.;

  Vector3i T;
  Vector3d v1;
  Vector3d v2;
  Vector3d v3;
  Vector3d gT;
  double s;

  for (int f=0 ; f<R.rows() ; f++){
    //we only add the triangles that belong to the region i
    if (R(f,0) == i){

      T = F.row(f);
      v1 = V.row(T(0));
      v2 = V.row(T(1));
      v3 = V.row(T(2));

      gT = g(v1,v2,v3);
      s = triangle_area(v1,v2,v3);

      Xi += s*gT;
      w += s;
    }
  }

  return Xi/w;

};

Vector3d new_Ni_L_2 (MatrixXi R, int i, MatrixXi F, MatrixXd V){

  //Compute the Covariance Matrix Ci
  MatrixXd Ci = MatrixXd::Zero(3,3); 
  double w;
  MatrixXd Xi = new_Xi_L_2(R,i,F,V);

  MatrixXd A(3,3);
  A(0,0) = 10; A(0,1) = 7; A(0,2) = 0;
  A(1,0) = 7; A(1,1) = 10; A(1,2) = 0;
  A(2,0) = 0; A(2,1) = 0.; A(2,2) = 0;

  Vector3i T;
  Vector3d v1;
  Vector3d v2;
  Vector3d v3;
  Vector3d gT;
  MatrixXd MT;
  double s;

  for (int f=0 ; f<R.rows() ; f++){
    //we only add the triangles that belong to the region i
    if (R(f,0) == i){

      T = F.row(f);
      v1 = V.row(T(0));
      v2 = V.row(T(1));
      v3 = V.row(T(2));

      gT = g(v1,v2,v3);
      MT = M(v1,v2,v3);
      s = triangle_area(v1,v2,v3);

      Ci += (2./72.)*s*MT*A*MT.transpose() + s*gT*gT.transpose();
      w += s;
    }
  }

  Ci = Ci - w*Xi*Xi.transpose();

  //Find the eigenvector of the min eigenvalue of Ci
  EigenSolver<MatrixXd> es(Ci);

  MatrixXd valp;
  MatrixXd vectp;
  valp = es.eigenvalues().real();
  vectp = es.eigenvectors().real();

  float min_valp;
  MatrixXd::Index minRow, minCol;
  Vector3d Ni;
  min_valp = valp.minCoeff(&minRow,&minCol);
  Ni = vectp.col(minRow);

  return Ni;

};

Vector3d new_Xi_L_2_1 (MatrixXi R, int i, MatrixXi F, MatrixXd V){

  return new_Xi_L_2(R,i,F,V);

};

Vector3d new_Ni_L_2_1 (MatrixXi R, int i, MatrixXi F, MatrixXd V){

  Vector3d Ni(0.,0.,0.);

  Vector3i T;
  Vector3d v1;
  Vector3d v2;
  Vector3d v3;
  double s;
  Vector3d nT;

  for (int f=0 ; f<R.rows() ; f++){
    //we only add the triangles that belong to the region i
    if (R(f,0) == i){

      T = F.row(f);
      v1 = V.row(T(0));
      v2 = V.row(T(1));
      v3 = V.row(T(2)); 

      s = triangle_area(v1,v2,v3);
      nT = triangle_normal(v1,v2,v3); 

      Ni += s*nT;
    }
  }

  return Ni;

};

//k is the number of regions/proxies of the partition R
MatrixXd new_proxies_L_2(MatrixXi R, MatrixXi F, MatrixXd V, int k){

  MatrixXd P(2*k,3);

  Vector3d Xi;
  Vector3d Ni;

  for (int i=0 ; i<k ; i++){
    Xi = new_Xi_L_2(R,i,F,V);
    Ni = new_Ni_L_2(R,i,F,V);
    P.row(i) = Xi;
    P.row(2*i) = Ni;
  }

  return P;

};

MatrixXd new_proxies_L_2_1(MatrixXi R, MatrixXi F, MatrixXd V, int k){

  MatrixXd P(2*k,3);

  Vector3d Xi;
  Vector3d Ni;

  for (int i=0 ; i<k ; i++){
    Xi = new_Xi_L_2_1(R,i,F,V);
    Ni = new_Ni_L_2_1(R,i,F,V);
    P.row(i) = Xi;
    P.row(2*i) = Ni;
  }

  return P;

};