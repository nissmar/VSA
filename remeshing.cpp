#include "remeshing.h"

bool list_contains(list <int> my_list, int my_var) {
    return (find(my_list.begin(), my_list.end(), my_var) != my_list.end());
}



int find_first_edge(HalfedgeDS he, int v, int r, MatrixXi R) {
    int edge = he.getOpposite(he.getEdge(v));
    while (R(he.getFace(edge),0)!=r || R(he.getFace(he.getOpposite(edge)),0)==r) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}
int find_next_edge(HalfedgeDS he, int edge, int r, MatrixXi R) {
    edge = he.getNext(edge);
    while (R(he.getFace(he.getOpposite(edge)),0)==r) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}

list <int> explore_boundary(HalfedgeDS he, int v, int r, MatrixXi R) {
    list <int> bound;
    int edge = find_first_edge(he,v,r,R);
    int new_v = he.getTarget(edge);
    bound.push_back(new_v);
    while (new_v != v) {
        edge = find_next_edge(he,edge,r,R);
        new_v = he.getTarget(edge);
        bound.push_back(new_v);
    }
    return bound;
}


int add_anchors_on_edge(HalfedgeDS he, int anchor, int r, MatrixXi R, MatrixXd V, VectorXi anchors) {
    // find next anchor
    int fedge = find_first_edge(he,anchor,r,R);
    int edge = fedge;
    int new_anchor = he.getTarget(edge);
    while (anchors(new_anchor)==0) {
        edge = find_next_edge(he,edge,r,R);
        new_anchor = he.getTarget(edge);
    }

    // find max between the two anchors
    cout << anchor<<" "<<new_anchor<<endl;
    edge = fedge;
    int new_v = he.getTarget(edge);
    double max_d = 0;
    int max_v = -1;
    while (new_v!= new_anchor) {
        edge = find_next_edge(he,edge,r,R);
        new_v = he.getTarget(edge);
        double dist = distance_projection(V,anchor, new_anchor,new_v);
        cout<<anchor<<" "<< new_anchor<<" "<<new_v<<" "<<dist<<endl;
        if (dist>max_d) {
            max_d = dist;
            max_v = new_v;
        }
    }
    return max_v;
}

list <int> find_vertex_proxies(HalfedgeDS he, int v, MatrixXi R) {
    list <int> vp;
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


VectorXi anchor_points(HalfedgeDS he, MatrixXi R, MatrixXd V) { // n is the number of vertices
    int n = V.rows();
    vector<list <int>> vertex_proxies(n); //list of proxies
    VectorXi anchors;
    anchors.setZero(n);
    int p=0; // number of anchors

    // find for each vertex its list of proxies
    for (int i=0;i<n;i++) {
        vertex_proxies[i] = find_vertex_proxies(he,i,R);
        if (vertex_proxies[i].size()>2) {
            anchors(i) = 1;
            p++; 

            // to show the boundary
            // int r =vertex_proxies[i].front();
            // list<int> b = explore_boundary(he,i,r,R);
            // cout << b.front() <<" oo"<<endl;
            // while(!b.empty()) {
            //     anchors(b.front()) = 1;
            //     b.pop_front();
            //     p++; 
            // }
        }
    }
    for (int i=0;i<n;i++) {
        vertex_proxies[i] = find_vertex_proxies(he,i,R);
        if (vertex_proxies[i].size()>2) {
            // to show the max deviation
            int r =vertex_proxies[i].front();

            int new_anchor= add_anchors_on_edge(he,i,r,R,V,anchors);
            cout << new_anchor<<endl;
            if (new_anchor!=-1) {
                anchors(new_anchor) = 1;
                p++;
            }
        }
    }
    VectorXi result;
    result.setZero(p); 
    for (int i=0;i<n;i++) {
        if (anchors(i)==1) {
            p--;
            result(p) = i;
        }
    }
    return result;
}

