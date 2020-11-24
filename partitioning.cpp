#include "partitioning.h"


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
  return fabs((N).dot(M-X))/N.norm();

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

double distance_L_2 (Vector3i T, Vector3d X, Vector3d N, MatrixXd V){

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

void tcolor(MatrixXi &Pf) {
  int n = Pf.rows();
  for (int i=0; i<n;i++) {
    Pf(i,0) = i;
  }
}

void fcolor(MatrixXd &Cf, MatrixXi Ad) {
  // closeness to the original face
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
  Cf.setZero(F.rows(),1);
  Vector3d c = triangle_center(F.row(0), V);
  Vector3d n = triangle_normal(F.row(0), V);
  for (int i=0;i<F.rows();i++) {
    Cf(i) = distance_L_2(F.row(i),c,n,V);
  }
  // closeness to the original face
}


void proxy_color(MatrixXi &Pf, MatrixXi Proxies, MatrixXd V, MatrixXi F, MatrixXi Ad) {
  int m=F.rows();
  Pf.setZero(m, 1);
  double color = 1.0;  


  priority_queue<pair<double, int>> q; // distance, face, proxy

  // initialize proxies
  int p = Proxies.rows();
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
    

    if (Pf(face)==0) {
      Pf(face) = prox;
      c++;
      // if (c>300) return;

      for (int k=0;k<3;k++) {
        int tri = Ad(face,k);
        double d = distance_L_2_1(F.row(tri), Proxies_normal[prox-1], V);
        q.push(make_pair(1.0-d, tri+m*prox));
      }
    }
  }
  
}

