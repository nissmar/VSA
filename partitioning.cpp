#include "partitioning.h"
#include "distance.h"
#include "proxies.h"

//******************Partitionning******************

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

MatrixXi uniform_proxies(int k, int n) {
  MatrixXi Proxies;
  Proxies.setZero(k,1);
  for (int i=0;i<k;i++) {
    Proxies(i) = (i*n)/k;
  }
  return Proxies;
}

void tcolor(MatrixXi &R) {
  // color according to the face number
  int n = R.rows();
  for (int i=0; i<n;i++) {
    R(i,0) = i;
  }
}

void fcolor(MatrixXd &Cf, MatrixXi Ad) {
  // closeness to the original face in terms of raw distance
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


void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V) {
  // closeness to the original face with the l2 or L21 distance
  Cf.setZero(F.rows(),1);
  Vector3d c = triangle_center(F.row(0), V);
  Vector3d n = triangle_normal(F.row(0), V);
  for (int i=0;i<F.rows();i++) {
    // Cf(i) = distance_L_2(F.row(i),c,n,V);
    Cf(i) = distance_L_2_1(F.row(i),n,V);
  }
}


void initial_partition(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad) {
  int m=F.rows();
  R.setZero(m, 1);

  priority_queue<pair<double, int>> q; // distance, face, proxy

  // initialize proxies
  MatrixXi Proxies = uniform_proxies(p, V.rows());
  Vector3d Proxies_center[p];
  Vector3d Proxies_normal[p];
  for (int i=0;i<p; i++) {
    Proxies_center[i] = triangle_center(F.row(Proxies(i)), V);
    Proxies_normal[i] = triangle_normal(F.row(Proxies(i)), V);
    q.push(make_pair(1.0,Proxies(i)+m*(i+1))); // proxies start at 1
  }

  pair<double, int> item;
  int face;
  int prox;

  int c=0;
  while (q.size()!=0) {
    item = q.top();
    q.pop();
    prox = item.second/m;
    face = item.second%m;
    

    if (R(face)==0) {
      R(face) = prox;
      c++;
      // if (c>300) return;

      for (int k=0;k<3;k++) {
        int tri = Ad(face,k);
        double d = distance_L_2(F.row(tri), Proxies_center[prox-1], Proxies_normal[prox-1], V);
        q.push(make_pair(1.0-d, tri+m*prox));
      }
    }
  }
  
}


VectorXi find_best_triangles(MatrixXd Proxies, MatrixXd V, MatrixXi F) {
  int p=Proxies.rows()/2;
  VectorXi triangles(p);
  VectorXd distances(p);
  double d;
  
  for (int i=0; i<p;i++) distances(i) = MAXFLOAT;

  for (int i=0; i<F.rows();i++) {
    for (int j=0; j<p;j++) {
      d = distance_L_2(F.row(i), Proxies.row(j), Proxies.row(j+p),  V);
      if (d < distances(j)) {
        distances(j) = d;
        triangles(j) = i;
      }
    }
  }
  return triangles;
}


void proxy_color(MatrixXi &R, MatrixXd Proxies, MatrixXd V, MatrixXi F, MatrixXi Ad) {
  int m=F.rows();
  R.setZero(m, 1);

  priority_queue<pair<double, int>> q; // distance, proxy

  // initialize proxies
  int p = Proxies.rows()/2;
  Vector3d Proxies_center[p];
  Vector3d Proxies_normal[p];
  VectorXi triangles = find_best_triangles(Proxies,V,F);
  for (int i=0;i<p; i++) {
    Proxies_center[i] = Proxies.row(i);
    Proxies_normal[i] = Proxies.row(i+p);
    q.push(make_pair(1.0, triangles(i)+m*(i+1))); // proxies start at 1
  }

  pair<double, int> item;
  int face;
  int prox;

  int c=0;
  while (q.size()!=0) {
    item = q.top();
    q.pop();
    prox = item.second/m;
    face = item.second%m;
    

    if (R(face)==0) {
      R(face) = prox;
      c++;
      for (int k=0;k<3;k++) {
        int tri = Ad(face,k);
        double d = distance_L_2(F.row(tri), Proxies_center[prox-1], Proxies_normal[prox-1], V);
        q.push(make_pair(1.0-d, tri+m*prox));
      }
    }
  }
  
}

