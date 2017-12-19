#include "dyy_RTK_method.hpp"
#include <stack>
#include <queue>

namespace dyy{

bool RTK::inkPW(dyy::RStarTree &tree, dyy::Mbr &Ew, int k, dyy::Point &q)
{
    std::queue<Node_P> queue;
    queue.push(tree.root);
    int rnk = 0;
    double score_qLow = dotMbrLow(q, Ew);
    double score_qUp = dotMbrUp(q, Ew);
    int rnkUp = 0;

    /*BFS*/
    while(!queue.empty()){
        Node_P e = queue.front();
        queue.pop();
        if(e->level){ //non leaf
            Node_P_V &children = *e->children;
            for(size_t ic = 0; ic < children.size(); ic++){
                Node_P childptr = children.at(ic);
                if(dotMMLow(childptr->mbrn, Ew) < score_qUp){
                    if(dotMMUp(childptr->mbrn, Ew) < score_qLow){
                        rnk += childptr->aggregate;
                        if(rnk >= k)
                            return false;
                    } else {
                        queue.push(childptr);
                    }
                }
            }
        } else { //leaf
            int rnk2 = rnk;
            Entry_P_V &entries = *e->entries;
            for(size_t ie = 0; ie < entries.size(); ie++){
                Entry_P entryPtr = entries.at(ie);
                double score_p_up = dotMMUp(entryPtr->mbre, Ew);
                //double score_p_low = dotMMLow(entryPtr->mbre, Ew);
                /*TBD dyy*/
                if(score_p_up < score_qUp){
                    if(score_p_up < score_qLow){
                        rnk++;
                        if(rnk >= k)
                            return false;
                    }
                    /*if(score_p_low < score_qUp){
                        rnkUp++;
                        }*/
                }
            }
        }
    }// BFS
    /*if(rnk + rnkUp < k)
        return true;
        else*/
    return false;
}

bool RTK::ink(dyy::RStarTree &tree, dyy::Point &w, int k, dyy::Point &q)
{
    std::queue<Node_P> queue;
    queue.push(tree.root);
    double score_q = dot(w,q);
    int rnk = 0;

    /*BFS tree*/
    while(!queue.empty()){
        Node_P e = queue.front();
        queue.pop();
        if(e->level){
            Node_P_V &children = *e->children;
            for(size_t ic = 0; ic < children.size(); ic++){
                Node_P childptr = children.at(ic);
                if(dotMbrLow(w,childptr->mbrn) < score_q){
                    if(dotMbrUp(w, childptr->mbrn) < score_q){
                        rnk += childptr->aggregate;
                        if(rnk >= k)
                            return false;
                    } else {
                        queue.push(childptr);
                    }
                }
            }
        } else { // leaf
            Entry_P_V &entries = *e->entries;
            for(size_t ie = 0; ie < entries.size();ie++){
                Entry_P entryPtr = entries.at(ie);
                //mbre is a point
                double score_p = dotMbrLow(w, entryPtr->mbre);
                if(score_p < score_q){
                    rnk++;
                    if(rnk >= k)
                        return false;
                }
            }
        }
    }

    return true;
}

Point_V RTK::rtkmethod(RStarTree &treeP, RStarTree &treeW, int k, Point &q)
{
    Point_V Results;

    std::stack<std::pair<Node_P, bool>> queue;

    queue.push(std::pair<Node_P, bool>(treeW.root, false));

    /*DFS Tree W*/
    while(!queue.empty()){
        Node_P e = queue.top().first;
        bool in = queue.top().second;
        queue.pop();

        if(e->level){ //non leaf
            if(in){
                Node_P_V &children = *e->children;
                for(size_t ic = 0; ic < children.size(); ic++){
                    Node_P childptr = children.at(ic);
                    queue.push(std::pair<Node_P, bool>(childptr,in));
                }
            } else {
                Node_P_V &children = *e->children;
                for(size_t ic = 0; ic < children.size(); ic++){
                    Node_P childptr = children.at(ic);
                    bool inFlag = RTK::inkPW(treeP, childptr->mbrn,k,q);
                    if(inFlag)
                        queue.push(std::pair<Node_P, bool>(childptr,true));
                    else
                        queue.push(std::pair<Node_P, bool>(childptr,false));
                }
            }
        } else { //leaf
            Entry_P_V &entries = *e->entries;
            for(size_t ie = 0; ie < entries.size(); ie++){
                Entry_P entryPtr = entries.at(ie);
                Point *w = static_cast<Point *>(entryPtr->data);
                if(in){
                    Results.push_back(*w);
                }else{
                    bool inFlag = RTK::ink(treeP, *w,k,q);
                    if(inFlag){
                        Results.push_back(*w);
                    }
                }
            }
        }
    }

    return Results;
}

}
