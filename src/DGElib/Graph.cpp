#include "Graph.h"
template<typename T>
//--------------------------------------------------------------------------
void Graph<T>::connexComp(typename std::vector<typename std::vector<T> > &grp) {
    std::vector<bool> mark(Graph<T>::vec_node.size(), false);

    for (typename std::vector<Graph<T>::Node *>::iterator i = Graph<T>::vec_node.begin();
         i != Graph<T>::vec_node.end(); i++) {
        if (!mark[(*i)->num_node])
            grp.push_back(bfs((*i)->num_node, mark));
    }
};

