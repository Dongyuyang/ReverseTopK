#ifndef _DYY_DATA_HPP_
#define _DYY_DATA_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include "const.h"
#include "dyy_RstarTree.hpp"

/*
  Data model used in the research of YuyangDong
 */

namespace dyy{
/*********************************************************************
 * Point
 ********************************************************************/
class Point
{
public:
    Coord_V coords;
    int id;
    Point(){};
    Point(Coord_V v){coords = v;}
    Point(Coord x, Coord y)
    {
        coords.push_back(x);
        coords.push_back(y);
    }

    friend std::istream& operator>>(std::istream& in, Point& p);
    bool operator< (const Point& p) const;

public:
    void print();
};


/*********************************************************************
 * Data
 ********************************************************************/

typedef std::vector<Point> Point_V;
typedef std::vector<LeafNodeEntry> Entry_V;

class Data
{
public:

    Point_V Products;
    Point_V Weights;
    Point_V Queries;

    Point_V Q2;

    RStarTree RtreeP;
    RStarTree RtreeW;

    Data(){};
    ~Data(){};

    Entry_V entriesP;
    Entry_V entriesW;

    static void loadPoint(std::string fileName, Point_V &v);
    static void buildTree(Point_V& points, Entry_V& entries, RStarTree* tree);

};


}



#endif /*_DYY_DATA_HPP_*/
