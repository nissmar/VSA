#include "proxies.h"
#include "distance.h"

Vector3d g(Vector3d v1,Vector3d v2,Vector3d v3){

  return (1./3.)*(v1+v2+v3);

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
  double w = 0.;
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
      
    }
    w += s;
  }
  
  Ci = Ci - w*Xi*Xi.transpose();
  
  //Find the eigenvector of the min eigenvalue of Ci
  EigenSolver<MatrixXd> es(Ci);

  MatrixXd valp;
  MatrixXd vectp;
  valp = es.eigenvalues().real();
  vectp = es.eigenvectors().real();

  double min_valp;
  MatrixXd::Index minRow, minCol;
  Vector3d Ni;
  min_valp = valp.minCoeff(&minRow,&minCol);
  Ni = vectp.col(minRow);

  return Ni.normalized();

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

  return Ni.normalized();

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
    P.row(k+i) = Ni.normalized();
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
    P.row(k+i) = Ni.normalized();
  }

  return P;

};

MatrixXd new_proxies(MatrixXi R, MatrixXi F, MatrixXd V, int k, int norme){
  if (norme == 0){
    return new_proxies_L_2(R,F,V,k);
  }
  else if (norme == 1){
    return new_proxies_L_2_1(R,F,V,k);
  }
  else {
    cout<<"wrong norme parameter"<<endl;
  }
};
