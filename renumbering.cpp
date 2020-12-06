#include "renumbering.h"

bool list_contains(list <int> my_list, int my_var) {
    return (find(my_list.begin(), my_list.end(), my_var) != my_list.end());
}

map<int,int> renumber(MatrixXi& F) {
    int n = F.rows();
    int m = F.cols();
    list<int> num;
    map<int,int> index;

    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) {
            int x = F(i,j);
            if (!list_contains(num,x)) {
                num.push_back(x);
            }
        }
    }
    num.sort();
    int i = 0;
    while (!num.empty()) {
        index[num.front()] = i;
        i++;
        num.pop_front();
    }
    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) {
            F(i,j) = index[F(i,j)];
        }
    }
    return index;
}

Vector3d projection(Vector3d z, MatrixXd Proxies, int k) {
    Vector3d p = Proxies.row(k);
    Vector3d n = Proxies.row(Proxies.rows()/2+k);
    double l = n.dot(p)-n.dot(z); //y is normalized
    return z-l*n;
}

MatrixXd new_V(HalfedgeDS he, MatrixXd V, MatrixXd Proxies, MatrixXi R, map<int,int> index){
    int ac = index.size(); // number of anchors
    MatrixXd newV;
    newV.setZero(ac,3);
    for (int i=0;i<V.rows();i++) {
        if (index.count(i)!=0) {
            vector<int> proxies = find_vertex_proxies(he,i,R);
            Vector3d p(0,0,0);
            for(size_t m = 0; m < proxies.size(); m++) {
                p += projection(V.row(i),Proxies,proxies[m]);
            }
            newV.row(index[i]) = p/proxies.size();
        }
    }
    return newV;
}