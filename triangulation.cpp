#include "triangulation.h"

int find_first(HalfedgeDS he, int v, int r, MatrixXi R) {
    int edge = he.getOpposite(he.getEdge(v));
    while (R(he.getFace(edge),0)!=r || R(he.getFace(he.getOpposite(edge)),0)==r) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}

int find_second(HalfedgeDS he, int v, int r, MatrixXi R) {
    int edge = he.getOpposite(he.getEdge(v));
    while (R(he.getFace(he.getOpposite(edge)),0)!=r || R(he.getFace(edge),0)==r) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}

int find_next_first(HalfedgeDS he, int edge, int r, MatrixXi R) {
    edge = he.getNext(edge);
    while (R(he.getFace(edge),0)!=r || R(he.getFace(he.getOpposite(edge)),0)==r) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}

int find_next_second(HalfedgeDS he, int edge, int r, MatrixXi R) {
    edge = he.getNext(edge);
    while (R(he.getFace(he.getOpposite(edge)),0)!=r || R(he.getFace(edge),0)==r) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}

vector<int> find_interior_neighbors(HalfedgeDS he, int v, int r, MatrixXi R){

    vector<int> neighbors;

    int edge = he.getEdge(v);
    //check if the vector is interior to the region r
    if (R(he.getFace(edge),0)==r || R(he.getFace(he.getOpposite(edge)),0)==r){
        neighbors.push_back(he.getOpposite(edge));
    }

    int edge_bis = he.getOpposite(he.getNext(edge));
    while (edge_bis != edge){
        if (R(he.getFace(edge_bis),0)==r || R(he.getFace(he.getOpposite(edge_bis)),0)==r){
            neighbors.push_back(he.getOpposite(edge_bis));
        }
        edge_bis = he.getOpposite(he.getNext(edge_bis));
    }

    return neighbors;

};

double get_length_edge(HalfedgeDS he, int edge, MatrixXd V){
    Vector3d v1 = V.row(he.getTarget(edge));
    Vector3d v2 = V.row(he.getTarget(he.getOpposite(edge)));
    return (v2-v1).norm();
};

MatrixXi color_region (MatrixXi R, int region, vector<vector<int>> anchors, MatrixXd V, HalfedgeDS he){

    priority_queue<pair<double, int>> q;
    MatrixXi color_graph = -MatrixXi::Ones(V.rows(),1); //each vertex of the region will be attributed a color, the others will stay at -1
    int nb_anchors = anchors[region].size();

    int anchor_vertex;
    int edge1;
    int edge2;
    //initiates a color for each anchor vertex of the region, and pushes in the queue its neighbors which are on the boundaries
    for (int i=0 ; i<nb_anchors ; i++){
        anchor_vertex = anchors[region][i];
        color_graph(anchor_vertex,0) = i;
        edge1 = find_first(he,anchor_vertex,region,R);
        edge2 = find_second(he,anchor_vertex,region,R);
        color_graph(he.getTarget(edge2)) = -2; //sign to recognise the direction 2
        q.push(make_pair(-get_length_edge(he,edge1,V),i+nb_anchors*edge1));
        q.push(make_pair(-get_length_edge(he,edge2,V),i+nb_anchors*edge2));
    }

    pair<double,int> item;
    int edge;
    int v;
    int color;
    double length;
    //empties the queue, assigns colors to boundary vertices and pushes their neighbors which are on the boundaries
    while(q.size()!=0){
        item = q.top();
        q.pop();
        edge = item.second/nb_anchors;
        v = he.getTarget(edge);
        color = item.second%nb_anchors;
        length = item.first;

        if (color_graph(v,0)==-1){
            color_graph(v,0) = color;
            edge = find_next_first(he,edge,region,R);
            q.push(make_pair(length-get_length_edge(he,edge,V),color+nb_anchors*edge));
        }
        
        else if (color_graph(v,0)==-2){
            color_graph(v,0) = color;
            edge = find_next_second(he,edge,region,R);
            q.push(make_pair(length-get_length_edge(he,edge,V),color+nb_anchors*edge));
        }
        
    }

    MatrixXi transition_color_graph = -MatrixXi::Ones(V.rows(),1); //for this phase, we need to go through boundary vertices so we ignore their color for now
    vector<int> edges;
    //does the same thing for all the interior vertices
    for (int i=0 ; i<nb_anchors ; i++){
        anchor_vertex = anchors[region][i];
        transition_color_graph(anchor_vertex,0) = i;
        edges = find_interior_neighbors(he,anchor_vertex,region,R);
        for (int j=0 ; j<edges.size() ; j++){
            q.push(make_pair(-get_length_edge(he,edges[j],V),i+nb_anchors*edges[j]));
        }
    }

    //empties the queue, assigns colors to interior vertices and pushes their neighbors which are inside the region
    while(q.size()!=0){
        item = q.top();
        q.pop();
        edge = item.second/nb_anchors;
        v = he.getTarget(edge);
        color = item.second%nb_anchors;
        length = item.first;

        if (transition_color_graph(v,0)==-1){
            transition_color_graph(v,0) = color;
            edges = find_interior_neighbors(he,he.getTarget(edge),region,R);
            for (int j=0 ; j<edges.size() ; j++){
                q.push(make_pair(length-get_length_edge(he,edges[j],V),color+nb_anchors*edges[j]));
            }
        }
        
    }

    //writing the interior vertices colors in the final graph
    for (int i=0 ; i<transition_color_graph.size() ; i++){
        if( (transition_color_graph(i,0) != -1) && (color_graph(i,0) == -1) ){
            color_graph(i,0) = transition_color_graph(i,0);
        }
    }

    /**for (int i=0 ; i<color_graph.size() ; i++){
        if(transition_color_graph(i,0) != -1){
            cout<<"vertex "<<i<<" color "<<color_graph(i,0)<<endl;
        }
    }*/

    return color_graph;

};

vector<Vector3i> triangulate_region (MatrixXi R, int region, vector<vector<int>> anchors, MatrixXd V, MatrixXi F, HalfedgeDS he){

    vector<Vector3i> triangles;

    MatrixXi color_graph = color_region(R,region,anchors,V,he);
    int nb_anchors = anchors[region].size();
    MatrixXi correspondence_color_anchor(nb_anchors,1);
    int color;
    int anchor_vertex;
    for (int i=0 ; i<nb_anchors ; i++){
        anchor_vertex = anchors[region][i];
        color = color_graph(anchor_vertex,0);
        correspondence_color_anchor(color,0) = anchor_vertex;
    }

    int color1;
    int color2;
    int color3;
    int vertex1;
    int vertex2;
    int vertex3;
    
    for (int f=0 ; f<R.size() ; f++){
        if (R(f,0)==region){
            color1 = color_graph(F(f,0));
            color2 = color_graph(F(f,1));
            color3 = color_graph(F(f,2));

            /**if (color1 == -1){
                cout<<"region "<<region<<" not fully colored, missing vertex "<<F(f,0)<<" face "<<f<<endl;
            }
            if (color2 == -1){
                cout<<"region "<<region<<" not fully colored, missing vertex "<<F(f,1)<<" face "<<f<<endl;
            }
            if (color3 == -1){
                cout<<"region "<<region<<" not fully colored, missing vertex "<<F(f,2)<<" face "<<f<<endl;
            }*/

            if (color1!=color2 && color1!=color3 && color2!=color3){
                vertex1 = correspondence_color_anchor(color1,0);
                vertex2 = correspondence_color_anchor(color2,0);
                vertex3 = correspondence_color_anchor(color3,0);
                triangles.push_back(Vector3i(vertex1,vertex2,vertex3));
            }
        }
    }

    return triangles;

};

pair<MatrixXi,MatrixXi> triangulation (MatrixXi R, vector<vector<int>> anchors, MatrixXd V, MatrixXi F, HalfedgeDS he){

    MatrixXi new_F;
    MatrixXi new_R;

    vector<Vector3i> triangulation;
    vector<int> regions;
    int nb_regions = anchors.size();

    vector<Vector3i> triangulation_region;
    for (int r=0 ; r<nb_regions ; r++){
        triangulation_region = triangulate_region(R,r,anchors,V,F,he);
        for (int i=0 ; i<triangulation_region.size() ; i++){
            triangulation.push_back(triangulation_region[i]);
            regions.push_back(r);
        }
    }

    int nb_triangles = triangulation.size();
    new_F = MatrixXi(nb_triangles,3);
    new_R = MatrixXi(nb_triangles,1);

    for (int i=0 ; i<nb_triangles ; i++){
        new_F.row(i) = triangulation[i];
        new_R(i,0) = regions[i];
    }

    return make_pair(new_F,new_R);

};