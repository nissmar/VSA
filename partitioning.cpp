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

vector<int> uniform_proxies(int k, int n) {
  vector<int> Proxies;
  for (int i=0;i<k;i++) {
    Proxies.push_back((i*n)/k);
  }
  return Proxies;
}


vector<int> random_proxies(int k, int n) {
  vector<int> Proxies;
  while (Proxies.size()<k){
    int x = rand() % n;  
    if (!vector_contains (Proxies,x)) Proxies.push_back(x);
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


void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V, int norme) {
  // closeness to the original face with the l2 or L21 distance
  Cf.setZero(F.rows(),1);
  Vector3d c = get_center(0);
  Vector3d n = get_normal(0);
  for (int i=0;i<F.rows();i++) {
    Cf(i) = distance(i,c,n,V, norme);
  }
}


int find_triangles_region(vector<int> Triangles, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad, int norme) {
  // updates R and find for each region the furthest triangle and its region
  int m=F.rows();
  R = -MatrixXi::Ones(m, 1);

  priority_queue<pair<double, int>> q; // distance, face, proxy
  

  // initialize proxies
  int p=Triangles.size();
  Vector3d Proxies_center[p];
  Vector3d Proxies_normal[p];

  // furthest 
  VectorXi furthest_triangle(p);
  VectorXd furthest_distance(p);
  furthest_distance.setZero(p);
  for (int i=0;i<p; i++) {
    Proxies_center[i] = get_center(Triangles[i]);
    Proxies_normal[i] = get_normal(Triangles[i]);
    R(Triangles[i]) = i;
    for (int k=0;k<3;k++) {
        int tri = Ad(Triangles[i],k);
        double d = distance(tri, Proxies_center[i], Proxies_normal[i], V, norme);
        q.push(make_pair(-d, tri+m*(i)));
        if (d>furthest_distance(i) && !vector_contains(Triangles,tri)) {
          furthest_distance(i)=d;
          furthest_triangle(i)=tri;
        }

    }
  }

  pair<double, int> item;
  int face;
  int prox;

  while (q.size()!=0) {
    item = q.top();
    q.pop();
    prox = item.second/m;
    face = item.second%m;
    

    if (R(face)==-1) {
      R(face) = prox;

      for (int k=0;k<3;k++) {
        int tri = Ad(face,k);
        double d = distance(tri, Proxies_center[prox], Proxies_normal[prox], V, norme);
        q.push(make_pair(-d, tri+m*prox));
        if (d>furthest_distance(prox) && !vector_contains(Triangles,tri)) {
          furthest_distance(prox)=d;
          furthest_triangle(prox)=tri;
        }
      }
    }
  }

  double max=furthest_distance(0);
  int maxtri=furthest_triangle(0);
  for (int i=1;i<p;i++) {
    if (max<furthest_distance(i)) {
      maxtri=furthest_triangle(i);
      max=furthest_distance(i);
    } 
  }
  return maxtri;
}

void initial_partition(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad, int norme) {
  vector<int> Triangles = random_proxies(p,F.rows());
  find_triangles_region(Triangles,R,V,F,Ad,norme);
}


void initial_partition2(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad, int norme) {
  //using furthest triangle
  int tri = rand() % F.rows();
  vector<int> Triangles;
  for (int i=0;i<p;i++) {
      cout << i<<endl;
      Triangles.push_back(tri);
      tri = find_triangles_region(Triangles,R,V,F,Ad,norme);
  }
}


VectorXi find_best_triangles(MatrixXi R, MatrixXd Proxies, MatrixXd V, MatrixXi F, int norme) {
  int p=Proxies.rows()/2;
  VectorXi triangles(p);
  VectorXd distances(p);
  double d;
  
  for (int i=0; i<p;i++) distances(i) = MAXFLOAT;
  
  //look at each face i
  for (int i=0; i<F.rows();i++) {
    int j = R(i,0);
    d = distance(i, Proxies.row(j), Proxies.row(j+p),  V, norme);
    if (d < distances(j)) {
      distances(j) = d;
      triangles(j) = i;
    }
  }
  return triangles;
}


double proxy_color(MatrixXi &R, MatrixXd Proxies, MatrixXd V, MatrixXi F, MatrixXi Ad, int norme) {
  int m=F.rows();
  double error=0;
  priority_queue<pair<double, int>> q; // distance, proxy

  // initialize proxies
  int p = Proxies.rows()/2;
  Vector3d Proxies_center[p];
  Vector3d Proxies_normal[p];
  VectorXi triangles = find_best_triangles(R,Proxies,V,F,norme);

  // reset R
  R = -MatrixXi::Ones(m, 1);
  for (int i=0;i<p; i++) {
    Proxies_center[i] = Proxies.row(i);
    Proxies_normal[i] = Proxies.row(i+p);
    R(triangles(i)) = i;
    error += distance(triangles(i), Proxies_center[i], Proxies_normal[i], V, norme);
    for (int k=0;k<3;k++) {
        int tri = Ad(triangles(i),k);
        double d = distance(tri, Proxies_center[i], Proxies_normal[i], V, norme);
        q.push(make_pair(-d, tri+m*(i)));
    }
  }

  pair<double, int> item;
  int face;
  int prox;

  while (q.size()!=0) {
    item = q.top();
    q.pop();
    error += - item.first;
    prox = item.second/m;
    face = item.second%m;
    

    if (R(face)==-1) {
      R(face) = prox;
      for (int k=0;k<3;k++) {
        int tri = Ad(face,k);
        double d = distance(tri, Proxies_center[prox], Proxies_normal[prox], V, norme);
        q.push(make_pair(-d, tri+m*prox));
      }
    }
  }
  return error;
}

