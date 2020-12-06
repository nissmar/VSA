#include "anchors.h"

bool list_contains(list <int> my_list, int my_var) {
    return (find(my_list.begin(), my_list.end(), my_var) != my_list.end());
}

bool vector_contains(vector<int> v, int key) {
    return count(v.begin(), v.end(), key);
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


int add_anchors_on_edge(HalfedgeDS he, int anchor, int r, MatrixXi R, MatrixXd V, VectorXi& anchors,  MatrixXd Proxies) {
    // modifies anchors to add anchors on edge
    int k = 0; //number of new anchors

    // find all anchors of the region R
    list<int>  anchors_list;
    int fedge = find_first_edge(he,anchor,r,R);
    int edge = fedge;
    int new_anchor = he.getTarget(edge);
    while (new_anchor != anchor) {
        if (anchors(new_anchor)==1) anchors_list.push_back(new_anchor);
        edge = find_next_edge(he,edge,r,R);
        new_anchor = he.getTarget(edge);
    }
    anchors_list.push_back(anchor);

    // go through each pair of anchor
    edge = fedge;
    int r2;
    while (!anchors_list.empty()) {
        new_anchor = anchors_list.front();
        anchors_list.pop_front();
        int new_v = he.getTarget(edge);
        double max_d = 0;
        int max_v = -1;
        while (new_v!= new_anchor) {
            edge = find_next_edge(he,edge,r,R);
            r2 = R(he.getFace(he.getOpposite(edge)),0);
            new_v = he.getTarget(edge);
            double dist = distance_projection(V, Proxies, anchor, new_anchor,new_v,r,r2);
            if (dist>max_d) {
                max_d = dist;
                max_v = new_v;
            }
        }
        if (max_d>0.4) { // treshold
            anchors(max_v) = 1;
            k++;
        }
        anchor = new_anchor;
    }
    return k;
}

vector<int> find_anchors_on_region(HalfedgeDS he, int anchor, int r, MatrixXi R, VectorXi anchors) {
    // modifies anchors to add anchors on edge
    int k = 0; //number of new anchors

    // find all anchors of the region R
    vector<int>  anchors_vector;
    anchors_vector.push_back(anchor);
    int fedge = find_first_edge(he,anchor,r,R);
    int edge = fedge;
    int new_anchor = he.getTarget(edge);
    while (new_anchor != anchor) {
        if (anchors(new_anchor)==1) {
            anchors_vector.push_back(new_anchor);
        }
        edge = find_next_edge(he,edge,r,R);
        new_anchor = he.getTarget(edge);
    }
    return anchors_vector;
}

vector<int> find_vertex_proxies(HalfedgeDS he, int v, MatrixXi R) {
    vector<int> vp;
    int proxy;
    int edge = he.getEdge(v);
    proxy = R(he.getFace(edge),0);
    vp.push_back(proxy);
    int pedge = he.getOpposite(he.getNext(edge));
    while (pedge != edge) {
        proxy = R(he.getFace(pedge),0);
        if (!vector_contains(vp,proxy)) {
            vp.push_back(proxy);
        }
        pedge = he.getOpposite(he.getNext(pedge));
    }
    return vp;
}


vector<vector<int>> anchor_points(HalfedgeDS he, MatrixXi R, MatrixXd V, MatrixXd Proxies) { 
    int n = V.rows();
    int p = Proxies.rows()/2;
    vector<vector<int>> vertex_proxies(n); //list of proxies
    VectorXi anchors;
    VectorXi seen; 
    seen.setZero(p);
    anchors.setZero(n);
    int k=0; // number of anchors

    // find for each vertex its list of proxies
    for (int i=0;i<n;i++) {
        vertex_proxies[i] = find_vertex_proxies(he,i,R);
        // add anchor if vertex has 3+ proxies
        if (vertex_proxies[i].size()>2) {
            anchors(i) = 1;
            k++; 
        }
    }

    // add anchors between existing ones
    int r,kv;
    for (int i=0;i<n;i++) {
        if (vertex_proxies[i].size()>2) {
            // to show the max deviation
            for(size_t m = 0; m < vertex_proxies[i].size(); m++) {
                r = vertex_proxies[i][m];
                if (seen(r)==0) {
                    seen(r)=1;
                    kv = add_anchors_on_edge(he,i,r,R,V,anchors, Proxies);
                    while (kv>0) {
                        k += kv;
                        kv = add_anchors_on_edge(he,i,r,R,V,anchors, Proxies);
                    }
                }
            }
        }
    }
    seen.setZero(p);
    // return list of polygons
    vector<vector<int>> polys_anchors;
    vector<int> void_vector;
    for (int i=0;i<p;i++) polys_anchors.push_back(void_vector);
    for (int i=0;i<n;i++) {
        if (vertex_proxies[i].size()>2) {
            for(size_t m = 0; m < vertex_proxies[i].size(); m++) {
                r = vertex_proxies[i][m];
                if (seen(r)==0) {
                    seen(r)=1;
                    polys_anchors[r] = find_anchors_on_region(he,i,r,R,anchors);
                }
            }
        }
    }
    return polys_anchors;
}

