#include "partitioning.h"


void tcolor(MatrixXi &Pf) {
  int n = Pf.rows();
  for (int i=0; i<n;i++) {
    Pf(i,0) = i;
  }
}
