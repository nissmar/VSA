#include "partitioning.h"

MatrixXi face_adjacency(MatrixXi F) {
  // face adjacency matrix
  int m = F.rows();
  MatrixXi Ad;
  Ad.setZero(m,3);
  int k;
  int corresp;
  for (int i=0;i<m;i++) {
    k = 0;
    for (int j=0;j<m;j++) {
      corresp = 0;
      for (int l=0; l<3; l++) {
        corresp += (F(i,l)==F(j,0)) + (F(i,l)==F(j,1)) + (F(i,l)==F(j,2));
      }
      if (corresp==2) {
        Ad(i,k) = j;
        k++;
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
  std::priority_queue<pair<double, int>> q;
  pair<double, int> face;
  q.push(make_pair(color,0));
  while (q.size()!=0) {
    face = q.top();
    q.pop();
    if (Cf(face.second,0)==0) {
      Cf(face.second,0) = face.first;
      color = face.first * 0.9;
      q.push(make_pair(color,Ad(face.second,0)));
      q.push(make_pair(color,Ad(face.second,1)));
      q.push(make_pair(color,Ad(face.second,2)));
    }
  }
  
}
