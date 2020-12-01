#include "remeshing.h"

bool triangle_contains(Vector3i T, int vertice) {
    return T(0)==vertice || T(1)==vertice || T(2)==vertice;
}

bool list_contains(list <int> my_list, int my_var) {
    return (find(my_list.begin(), my_list.end(), my_var) != my_list.end());
}

bool is_anchor(MatrixXi R, MatrixXi Ad, MatrixXi F, int face, int vertice) {
    list <int> seen, to_see;
    to_see.push_back(face);
    int tri;
    int first_color = R(face,0);
    int second_color = -1;
    int color;
    while (!to_see.empty()) {
        face = to_see.back();
        to_see.pop_back();
        if (!list_contains(seen,face)){
            seen.push_back(face);
            for (int k=0;k<3;k++) {
                tri = Ad(face,k);
                color = R(face,0);
                if (triangle_contains(F.row(tri),vertice)) {
                    if (color!=first_color) {
                        if (second_color==-1) second_color = color;
                        if (color!=second_color) return true;
                    }
                    to_see.push_back(tri);
                }
            }
        }
    }
    return false;
}

VectorXi anchor_points(MatrixXi R,  MatrixXi Ad, MatrixXi F, int n) { // n is the number of vertices
    VectorXi anchors;
    int p=0;
    anchors.setZero(n); // 0 = not seen 1 = not anchor 2 = anchor
    for (int i=0;i<F.rows();i++) {
        for (int k=0;k<3;k++) {
            int vertice = F(i,k);
            if (anchors(vertice)==0) {
                if (is_anchor(R, Ad,  F, i, vertice)) {
                    anchors(vertice)=2;
                    p++;
                }
                else anchors(vertice)=1;
            }
        }
    }
    VectorXi result;
    result.setZero(p); // 0 = not seen 1 = not anchor 2 = anchor
    for (int i=0;i<n;i++) {
        if (anchors(i)==2) {
            p--;
            result(p) = i;
        }
    }
    return result;

}