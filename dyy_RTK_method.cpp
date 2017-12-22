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


RTK::ARankResult RTK::inRank(RStarTree &tree, Point &w, int minRank, Point &q)
{
    std::queue<Node_P> queue;
    queue.push(tree.root);
    int rnk = 0;
    double q_score = dot(q, w);
    ARankResult AR;
    AR.isBetter = false;

    /*BFS*/
    while(!queue.empty()){
        Node_P e = queue.front();
        queue.pop();
        if(e->level){ //non leaf
            Node_P_V &children = *e->children;
            for(size_t ic = 0; ic < children.size(); ic++){
                Node_P childptr = children.at(ic);
                if(dotMbrLow(w,childptr->mbrn) < q_score){
                    if(dotMbrUp(w, childptr->mbrn) < q_score){
                        rnk += childptr->aggregate;
                        if(rnk >= minRank)
                            return AR;
                    } else {
                        queue.push(childptr);
                    }
                }
            }
        } else { //leaf
            Entry_P_V &entries = *e->entries;
            for(size_t ie = 0; ie < entries.size();ie++){
                Entry_P entryPtr = entries.at(ie);
                //mbre is a point
                double score_p = dotMbrLow(w, entryPtr->mbre);
                if(score_p < q_score){
                    rnk++;
                    if(rnk >= minRank)
                        return AR;
                }
            }
        }
    }// while BFS search

    if(rnk < minRank){
        AR.isBetter = true;
        AR.rank = rnk;
    }

    return AR;
}



RTK::ARankResult RTK::inRankPW(RStarTree &tree, Mbr &Ew, int minRank, Point &q)
{
        std::queue<Node_P> queue;
    queue.push(tree.root);
    int rnk = 0;
    double score_qLow = dotMbrLow(q, Ew);
    double score_qUp = dotMbrUp(q, Ew);

    ARankResult AR;
    AR.flag = -1;

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
                        if(rnk >= minRank)
                            return AR;
                    } else {
                        queue.push(childptr);
                    }
                }
            }
        } else { //leaf
            Entry_P_V &entries = *e->entries;
            for(size_t ie = 0; ie < entries.size(); ie++){
                Entry_P entryPtr = entries.at(ie);
                double score_p_up = dotMMUp(entryPtr->mbre, Ew);
                if(score_p_up < score_qUp){
                    if(score_p_up < score_qLow){
                        rnk++;
                        if(rnk >= minRank)
                            return AR;
                    }
                }
            }
        }
    }// BFS

    /*not sure*/
    if(rnk < minRank)
        AR.flag = 1;
    else
        AR.flag = 0;
    return AR;
}


RTK::BUFFER RTK::rkrmethod(RStarTree &treeP, RStarTree &treeW, Point &q, int k)
{
    BUFFER buffer;
    int threshold = std::numeric_limits<int>::max();
    std::stack<std::pair<Node_P, bool>> queue;

    queue.push(std::pair<Node_P,bool>(treeW.root,false));


    /*DFS W*/
    while(!queue.empty()){
        Node_P e = queue.top().first;
        bool in = queue.top().second;
        queue.pop();
        if(buffer.size() == k)
            threshold = std::prev(buffer.end())->first;
        if(e->level){ //non leaf
            if(in){
                Node_P_V &children = *e->children;
                for(size_t ic = 0; ic < children.size(); ic++){
                    Node_P childptr = children.at(ic);
                    queue.push(std::pair<Node_P, bool>(childptr,in));
                }

            }
            else{
            Node_P_V &children = *e->children;
            for(size_t ic = 0; ic < children.size(); ic++){
                Node_P childptr = children.at(ic);
                if(threshold == std::numeric_limits<int>::max()){
                    queue.push(std::pair<Node_P, bool>(childptr,false));
                }
                else{
                    ARankResult ar = inRankPW(treeP, childptr->mbrn, threshold, q);
                    if(ar.flag == 0)
                        queue.push(std::pair<Node_P, bool>(childptr,false));
                    if(ar.flag == 1)
                        queue.push(std::pair<Node_P, bool>(childptr,true));
                }
            }
            }
        } else { //leaf
            Entry_P_V &entries = *e->entries;
            for(size_t ie = 0; ie < entries.size(); ie++){
                Entry_P entryPtr = entries.at(ie);
                Point *w = static_cast<Point *>(entryPtr->data);
                ARankResult ar = inRank(treeP, *w, threshold, q);
                if(ar.isBetter){
                    if(buffer.size() == k)
                        update_map(buffer, ar.rank, w->id);
                    else
                        buffer.insert(std::pair<int, int>(ar.rank, w->id));
                }
            }
        }
    }// loop BFS W

    return buffer;
}

}
