#include "remeshing.h"

bool list_contains(list <int> my_list, int my_var) {
    return (find(my_list.begin(), my_list.end(), my_var) != my_list.end());
}

list <int> vertex_proxies(HalfedgeDS he, int v, MatrixXi R) {
    list <int> vp;
    // TO BE COMPLETED
    int proxy;
    int edge = he.getEdge(v);
    proxy = R(he.getFace(edge),0);
    vp.push_back(proxy);
    int pedge = he.getOpposite(he.getNext(edge));
    while (pedge != edge) {
        proxy = R(he.getFace(pedge),0);
        if (!list_contains(vp,proxy)) {
            vp.push_back(proxy);
        }
        pedge = he.getOpposite(he.getNext(pedge));
    }
    return vp;
}


VectorXi anchor_points(HalfedgeDS he, MatrixXi R, int n) { // n is the number of vertices
    VectorXi anchors;
    int p=0;
    anchors.setZero(n); // 0 = not seen 1 = not anchor 2 = anchor
    list <int> vp;
    for (int i=0;i<n;i++) {
        vp = vertex_proxies(he,i,R);
        if (vp.size()>2) {
            anchors(i)=2;
            p++; 
        }
        else anchors(i)=1;
    }
    VectorXi result;
    result.setZero(p); 
    for (int i=0;i<n;i++) {
        if (anchors(i)==2) {
            p--;
            result(p) = i;
        }
    }
    return result;
}

