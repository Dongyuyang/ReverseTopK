#ifndef _DYY_RSTARTREE_HPP_
#define _DYY_RSTARTREE_HPP_

#include <algorithm>
#include <list>
#include "const.h"

#define MAX_CHILD (pageSize - sizeof(RTreeNode)) / sizeof(RTreeNode)


namespace dyy
{

const size_t LEAF_LEVEL = 0;
const size_t PAGE_SIZE = 4096;

const double SPLIT_FACTOR = 0.4;
const double REINSERT_FACTOR = 0.3;
const size_t NEAR_MINIMUM_OVERLAP_FACTOR = 32;

/******************************************************************************
 * MBR
 *****************************************************************************/

class Mbr
{
public:
    Coord_VV coord;
    Mbr();
    Mbr(Coord_V v);
    Mbr(Coord_VV vv);
    Mbr(Coord xmin, Coord xmax, Coord ymin, Coord ymax);

    void init();
    void print();
    double getArea() const;
    double getMargin() const;
    void enlarge(Mbr &add);

    inline double getCenter(size_t dim)
    {return (coord[dim][0] + coord[dim][1]) * 0.5;}

    static Mbr getMbr(Mbr &mbr1, Mbr &mbr2);
    static double getOverlap(Mbr &mbr1, Mbr &mbr2);
    static double getEnlarge(Mbr &base, Mbr &add);
};

/******************************************************************************
 * RtreeNode (leaf, noneleaf)
 *****************************************************************************/

typedef void *Data_P;
/*Leaf*/
class LeafNodeEntry
{
public:
    Mbr mbre;
    Data_P data;
    double value;

    LeafNodeEntry(Mbr _mbre, Data_P _data) : mbre(_mbre), data(_data){}
    void print(std::string ident);
};

/*NonLeaf*/
class RTreeNode;/*cliam class*/
typedef LeafNodeEntry* Entry_P;
typedef std::vector<Entry_P> Entry_P_V;
typedef RTreeNode* Node_P;
typedef std::vector<Node_P> Node_P_V;
typedef std::list<Node_P> Node_P_L;


class RTreeNode
{
public:
    size_t level;// level of node, leaf node is 0;
    size_t aggregate;
    Mbr mbrn;
    Node_P_V* children;
    Entry_P_V* entries;
    Node_P parent;
    double value;

    RTreeNode(size_t _level);
    ~RTreeNode();

    bool operator < (const RTreeNode& node) const {return value < node.value;}

    void print(std::string ident);

    size_t size(){ return level ? children->size() : entries->size();}

    void insert(Node_P childrenPtr);
    void insert(Entry_P entryPtr);

    void take(RTreeNode &source, size_t from);
    void resize(size_t from);
    void adjustMbr();

    void enlargeMbr(Mbr &mbr);
    void decreaseAggregate(size_t amount);
    void increaseAggregate(size_t amount);

    void prepareEnlargementSort(Mbr &add);
    void prepareDisSort();
};



/******************************************************************************
 * Sorting
 *****************************************************************************/

class AxisSort
{
public:
    size_t axis;

    AxisSort(size_t _axis) : axis(_axis) {};

    bool operator() (Node_P child1, Node_P child2);
    bool operator() (Entry_P entry1, Entry_P entry2);
};




/******************************************************************************
 * RStarTree
 *****************************************************************************/

class RStarTree{
public:
    static size_t pageSize; //Page size in bytes.
    static size_t maxChild; //maximum number of entries in a node.
    static size_t minChild; //minimum number of entries in a node.

    RTreeNode* root;

    RStarTree() : root(NULL) { root = new RTreeNode(LEAF_LEVEL);}
    ~RStarTree() {
        if(root)
            delete root;
    };

    void print();

    /// Insertion
    void insertData(Entry_P entryPtr);

private:
    /// Insertion
    /// with OverflowTreatment
    void insert(Node_P childPtr, Entry_P entryPtr,
                size_t desiredLevel, size_t &overflowLevel);

    /// Insertion - ChooseSubtree
    Node_P chooseSubtree(Mbr &mbr, size_t desiredLevel);

    /// Insertion - OverflowTreatment - ReInsert
    void reInsert(RTreeNode &node, size_t &overflowLevel);

    /// Insertion - OverflowTreatment - Split
    void split(RTreeNode &node);
    size_t chooseSplitAxis(RTreeNode &node);
    double computeS(RTreeNode &node);
    size_t chooseSplitIndex(RTreeNode &node, size_t axis);
};




}


#endif /*_DYY_RSTARTREE_HPP_*/
