#include "dyy_RstarTree.hpp"
#include <iostream>

namespace dyy
{

/******************************************************************************
 * MBR
 *****************************************************************************/

Mbr::Mbr()
{
    init();
}

void Mbr::init()
{
    Coord_V v = {INF_P, INF_N};
    for(size_t dim = 0; dim < DIM; dim++)
        coord.push_back(v);
}

Mbr::Mbr(Coord_V v)
{
    for(size_t dim = 0; dim < DIM; dim++)
        coord.push_back({v[dim], v[dim]});
}

Mbr::Mbr(Coord_VV vv)
{
    coord = vv;
}

Mbr::Mbr(Coord xmin, Coord xmax, Coord ymin, Coord ymax)
{
    coord.push_back({xmin,xmax});
    coord.push_back({ymin,ymax});
}

void Mbr::print()
{
    for(size_t dim = 0; dim < DIM; dim++)
        std::cout << "[" << coord[dim][0] << "," << coord[dim][1] << "]";
    std::cout << std::endl;
}

Mbr Mbr::getMbr(Mbr &mbr1, Mbr &mbr2)
{
    Mbr mbr;
    for(size_t dim = 0; dim < DIM; dim++){
        mbr.coord[dim][0] = std::min(mbr1.coord[dim][0], mbr2.coord[dim][0]);
        mbr.coord[dim][1] = std::max(mbr1.coord[dim][1], mbr2.coord[dim][1]);
    }
    return mbr;
}

double Mbr::getOverlap(Mbr &mbr1, Mbr &mbr2)
{
    double overlap = 1;
    for(size_t dim = 0; dim < DIM; dim++){
        Coord maxMin = std::max(mbr1.coord[dim][0], mbr2.coord[dim][0]);
        Coord minMax = std::min(mbr1.coord[dim][1], mbr2.coord[dim][1]);
        if(maxMin >= minMax)
            return 0;
        overlap *= minMax - maxMin;
    }
    return overlap;
}

double Mbr::getEnlarge(Mbr &base, Mbr &add)
{
    Mbr mbr = Mbr::getMbr(base, add);
    return mbr.getArea() - base.getArea();
}

double Mbr::getArea() const
{
    double area = 1;
    for(size_t dim = 0; dim < DIM; dim++)
        area *= coord[dim][1] - coord[dim][0];
    return area;
}

double Mbr::getMargin() const
{
    double margin = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        margin += coord[dim][1] - coord[dim][0];
    return margin;
}

void Mbr::enlarge(Mbr &add)
{
    for(size_t dim = 0; dim < DIM; dim++){
        coord[dim][0] = std::min(coord[dim][0], add.coord[dim][0]);
        coord[dim][1] = std::max(coord[dim][1], add.coord[dim][1]);
    }
}


/******************************************************************************
 * RtreeNode (leaf, noneleaf)
 *****************************************************************************/


void LeafNodeEntry::print(std::string ident)
{
    std::cout << ident;
    mbre.print();
}

RTreeNode::RTreeNode(size_t _level)
    : level(_level), aggregate(0), children(NULL), entries(NULL), parent(NULL)
{
    if(level) // if nonleaf
        children = new Node_P_V();
    else
        entries = new Entry_P_V();
}

RTreeNode::~RTreeNode()
{
    if(children){
        for(size_t ic = 0; ic < children->size(); ic++)
            delete children->at(ic);
        delete children;
    }
    if(entries)
        delete entries;
}

void RTreeNode::print(std::string ident)
{
    std::cout << ident << level << " " << aggregate;
    mbrn.print();
    if(level){
        for(size_t ic = 0; ic < children->size(); ic++)
            children->at(ic)->print(ident + "| ");
    } else {
        for(size_t ie = 0; ie < entries->size(); ie++)
            entries->at(ie)->print(ident + "| ");
    }
}

void RTreeNode::insert(Node_P childrenPtr)
{
    children->push_back(childrenPtr);
    childrenPtr->parent = this;
    enlargeMbr(childrenPtr->mbrn);
    increaseAggregate(childrenPtr->aggregate);

}


void RTreeNode::insert(Entry_P entryPtr)
{
    entries->push_back(entryPtr);
    enlargeMbr(entryPtr->mbre);
    increaseAggregate(1);
}


void RTreeNode::take(RTreeNode &source, size_t from)
{
    if(source.level){ //nonleaf
        Node_P_V &sourceV = *source.children;
        for(size_t ic = from; ic < sourceV.size(); ic++)
            insert(sourceV.at(ic));
    } else { //leaf
        Entry_P_V &sourceV = *source.entries;
        for(size_t ic = from; ic < sourceV.size(); ic++)
            insert(sourceV.at(ic));
    }
    source.resize(from);
}

void RTreeNode::resize(size_t from)
{
    if(level){
        children->resize(from);
        size_t newAggregate = 0;
        for(size_t ic = 0; ic < children->size(); ic++)
            newAggregate += children->at(ic)->aggregate;
        decreaseAggregate(aggregate - newAggregate);
    } else {
        entries->resize(from);
        decreaseAggregate(aggregate - from);
    }
    adjustMbr();
}

void RTreeNode::adjustMbr()
{
    mbrn.init();
    if(level)
        for(size_t ic = 0; ic < children->size(); ic++)
            mbrn.enlarge(children->at(ic)->mbrn);
    else
        for(size_t ie = 0; ie < entries->size(); ie++)
            mbrn.enlarge(entries->at(ie)->mbre);

    if(parent)
        parent->adjustMbr();
}


void RTreeNode::enlargeMbr(Mbr &mbr)
{
    mbrn.enlarge(mbr);
    if(parent)
        parent->enlargeMbr(mbr);
}

void RTreeNode::decreaseAggregate(size_t amount)
{
    aggregate -= amount;
    if(parent)
        parent->decreaseAggregate(amount);
}

void RTreeNode::increaseAggregate(size_t amount)
{
    aggregate += amount;
    if(parent)
        parent->increaseAggregate(amount);
}

void RTreeNode::prepareEnlargementSort(Mbr &add)
{
    for(size_t ic = 0; ic < children->size(); ic++){
        Node_P childPtr = children->at(ic);
        childPtr->value = Mbr::getEnlarge(childPtr->mbrn, add);
    }
}

void RTreeNode::prepareDisSort()
{
    Coord center[DIM];
    for(size_t dim = 0; dim < DIM; dim++)
        center[dim] = mbrn.getCenter(dim);
    if(level){
        for(size_t ic = 0; ic < children->size(); ic++){
            Node_P childPtr = children->at(ic);
            childPtr->value = 0;
            for(size_t dim = 0; dim < DIM; dim++){
                double sideLength = childPtr->mbrn.getCenter(dim) - center[dim];
                childPtr->value += sideLength * sideLength;
            }
        }
    } else {
        for(size_t ie = 0; ie < entries->size(); ie++){
            Entry_P entryPtr = entries->at(ie);
            entryPtr->value = 0;
            for(size_t dim = 0; dim < DIM; dim++){
                double sideLength = entryPtr->mbre.getCenter(dim) - center[dim];
                entryPtr->value += sideLength * sideLength;
            }
        }
    }
}


/******************************************************************************
 * Sorting
 *****************************************************************************/

bool AxisSort::operator()(Node_P child1, Node_P child2)
{
    if(child1->mbrn.coord[axis][0] == child2->mbrn.coord[axis][0])
        return child1->mbrn.coord[axis][1] < child2->mbrn.coord[axis][1];
    return child1->mbrn.coord[axis][0] < child2->mbrn.coord[axis][0];
}

bool AxisSort::operator()(Entry_P entry1, Entry_P entry2)
{
    if(entry1->mbre.coord[axis][0] == entry2->mbre.coord[axis][0])
        return entry1->mbre.coord[axis][1] < entry2->mbre.coord[axis][1];
    return entry1->mbre.coord[axis][0] < entry2->mbre.coord[axis][0];
}


/******************************************************************************
 * RStarTree
 *****************************************************************************/
size_t RStarTree::pageSize = PAGE_SIZE;
size_t RStarTree::maxChild = MAX_CHILD;
size_t RStarTree::minChild = maxChild * SPLIT_FACTOR;


void RStarTree::print()
{
    root->print("");
}


/// Invoke Insert() starting with the leaf level as a parameter,
/// to insert a new data rectangle
void RStarTree::insertData(Entry_P entryPtr)
{
    size_t overflowLevel = -1;
    overflowLevel <<= (root->level);
    insert(NULL, entryPtr, LEAF_LEVEL, overflowLevel);
}


void RStarTree::insert(Node_P childPtr, Entry_P entryPtr,
                       size_t desiredLevel, size_t &overflowLevel)

{
    Mbr &mbr = (desiredLevel ? childPtr->mbrn : entryPtr->mbre);
    Node_P node = chooseSubtree(mbr, desiredLevel);
    if(desiredLevel)
        node->insert(childPtr);
    else
        node->insert(entryPtr);

    while(node){ // while node is not NULL
        while(node->size() > RStarTree::maxChild){
            if(overflowLevel & (1 << node->level))
                split(*node);
            else
                reInsert(*node, overflowLevel);
        }
        node = node->parent;
    }
}


Node_P RStarTree::chooseSubtree(Mbr &mbr, size_t desiredLevel)
{
    Node_P node = root;
    while(node->level > desiredLevel){
        /// since (node->level > desiredlevel), node->level > 0, node is not leaf
        Node_P_V &children = *node->children;
        //node->
        node->prepareEnlargementSort(mbr);
        std::sort(children.begin(), children.end());

        size_t selectedIndex = 0;
        if(node->level == 1){ // if the child pointers point to leaf nodes
            size_t p = std::min(node->size(), NEAR_MINIMUM_OVERLAP_FACTOR);
            double minOverlapEnlarge = INF_P;
            for(size_t ic = 0; ic < p; ic++){
                Mbr base = children.at(ic)->mbrn;
                Mbr newMbr = Mbr::getMbr(base, mbr);
                double overlapBase = 0, overlap = 0;
                for(size_t ico = 0; ico < node->size(); ico++){
                    if(ico != ic){
                        overlapBase += Mbr::getOverlap(base, children.at(ico)->mbrn);
                        overlap += Mbr::getOverlap(newMbr, children.at(ico)->mbrn);
                    }
                }
                double overlapEnlarge = overlap - overlapBase;
                if(overlapEnlarge - EPS < minOverlapEnlarge){
                    minOverlapEnlarge = overlapEnlarge;
                    selectedIndex = ic;
                }
            }
        }

        node = children.at(selectedIndex);
    }

    return node;
}

void RStarTree::reInsert(RTreeNode &node, size_t &overflowLevel)
{
    overflowLevel |= (1 << node.level);

    node.prepareDisSort();
    if(node.level)
        std::sort(node.children->begin(), node.children->end());
    else
        std::sort(node.entries->begin(), node.entries->end());

    size_t p = RStarTree::maxChild * REINSERT_FACTOR;
    size_t reinsertFromIndex = node.size() - p;

    if(node.level){
        Node_P_V children;
        for(size_t ic = reinsertFromIndex; ic < node.children->size(); ic++)
            children.push_back(node.children->at(ic));
        node.resize(reinsertFromIndex);
        for(size_t ic = 0; ic < children.size();ic++)
            insert(children.at(ic), NULL, node.level, overflowLevel);
    } else {
        Entry_P_V entries;
        for(size_t ie = reinsertFromIndex; ie < node.entries->size(); ie++)
            entries.push_back(node.entries->at(ie));
        node.resize(reinsertFromIndex);
        for(size_t ie = 0; ie < entries.size(); ie++)
            insert(NULL, entries.at(ie), node.level, overflowLevel);
    }
}


void RStarTree::split(RTreeNode &node)
{
    size_t axis = chooseSplitAxis(node);
    size_t splitIndex = chooseSplitIndex(node, axis);

    Node_P newNode = new RTreeNode(node.level);
    newNode->take(node,splitIndex);
    if(!node.parent){
        root = new RTreeNode(node.level + 1);
        root->insert(&node);
    }
    node.parent->insert(newNode);
}


size_t RStarTree::chooseSplitAxis(RTreeNode &node)
{
    double minS = INF_P;
    size_t axis = 0;
    for(size_t dim = 0; dim < DIM; dim++){
        if(node.level)
            std::sort(node.children->begin(), node.children->end(), AxisSort(dim));
        else
            std::sort(node.entries->begin(), node.entries->end(), AxisSort(dim));
        double S = computeS(node);
        if(S - EPS < minS){
            minS = S;
            axis = dim;
        }
    }
    return axis;
}

double RStarTree::computeS(RTreeNode &node)
{
    double S = 0;
    size_t size = node.size(), last = size - 1;

    std::vector<Mbr> partition1(size), partition2(size);
    //Mbr partition1[size], partition2[size];
    if(node.level){
        Node_P_V &children = *node.children;
        partition1[0] = Mbr::getMbr(children.at(0)->mbrn, children.at(0)->mbrn);
        for(size_t ic = 1; ic < size - minChild; ic++)
            partition1[ic] = Mbr::getMbr(partition1[ic - 1], children.at(ic)->mbrn);
        partition2[last] = Mbr::getMbr(children.at(last)->mbrn, children.at(last)->mbrn);
        for(size_t ic = last - 1; ic >= minChild; ic--)
            partition2[ic] = Mbr::getMbr(partition2[ic + 1], children.at(ic)->mbrn);

    } else {
        Entry_P_V &entries = *node.entries;
        partition1[0] = Mbr::getMbr(entries.at(0)->mbre, entries.at(0)->mbre);
        for(size_t ie = 1; ie < size - minChild; ie++)
            partition1[ie] = Mbr::getMbr(partition1[ie - 1], entries.at(ie)->mbre);
        partition2[last] = Mbr::getMbr(entries.at(last)->mbre, entries.at(last)->mbre);
        for(size_t ie = last - 1; ie >= minChild; ie--)
            partition2[ie] = Mbr::getMbr(partition2[ie + 1], entries.at(ie)->mbre);
    }
    // is : first element of a valid second partition
    for(size_t is = minChild; is <= size - minChild; is++)
        S += partition1[is - 1].getMargin() + partition2[is].getMargin();

    return S;
}



size_t RStarTree::chooseSplitIndex(RTreeNode &node, size_t axis)
{
    size_t size = node.size(), last = size - 1;

    std::vector<Mbr> partition1(size), partition2(size);
    if(node.level){
        std::sort(node.children->begin(), node.children->end(), AxisSort(axis));
        Node_P_V &children = *node.children;
        partition1[0] = Mbr::getMbr(children.at(0)->mbrn, children.at(0)->mbrn);
        for(size_t ic = 1; ic < size - minChild; ic++)
            partition1[ic] = Mbr::getMbr(partition1[ic - 1], children.at(ic)->mbrn);
        partition2[last] = Mbr::getMbr(children.at(last)->mbrn, children.at(last)->mbrn);
        for(size_t ic = last - 1; ic >= minChild; ic--)
            partition2[ic] = Mbr::getMbr(partition2[ic + 1], children.at(ic)->mbrn);
    } else {
        std::sort(node.entries->begin(), node.entries->end(), AxisSort(axis));
        Entry_P_V &entries = * node.entries;
        partition1[0] = Mbr::getMbr(entries.at(0)->mbre, entries.at(0)->mbre);
        for(size_t ie = 1; ie < size - minChild; ie++)
            partition1[ie] = Mbr::getMbr(partition1[ie - 1], entries.at(ie)->mbre);
        for(size_t ie = last - 1; ie >= minChild; ie--)
            partition2[ie] = Mbr::getMbr(partition2[ie + 1], entries.at(ie)->mbre);
    }

    double minOverlap = INF_P, minArea = INF_P;
    size_t splitIndex = minChild;
    //is: first elemnt of a vaild second partition
    for(size_t is = minChild; is <= size - minChild; is++){
        double overlap = Mbr::getOverlap(partition1[is - 1], partition2[is]);
        if(overlap - EPS < minOverlap){
            double area = partition1[is - 1].getArea() + partition2[is].getArea();
            if(overlap + EPS > minOverlap){
                if(area < minArea){
                    minArea = area;
                    splitIndex = is;
                }
            } else {
                minOverlap = overlap;
                minArea = area;
                splitIndex = is;
            }
        }
    }

    return splitIndex;
}






}/*namespace: dyy*/

