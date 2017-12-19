#ifndef _DYY_RTK_METHOD_HPP_
#define _DYY_RTK_METHOD_HPP_

#include "dyy_data.hpp"
#include <map>
#include <vector>
#include <iterator>

namespace dyy{

class IntopK
{
public:
    bool intopk;
};


class RTK
{
public:

    /*dot procut*/
    static double dot(Point &a, Point &b);
    static double dotMbrLow(Point &point, Mbr &mbr);
    static double dotMbrUp(Point &point, Mbr &mbr);
    static double dotMMLow(Mbr &a, Mbr &b);
    static double dotMMUp(Mbr &a, Mbr &b);


    static bool inkPW(RStarTree &tree, Mbr &Ew, int k, Point &q);
    static bool ink(RStarTree &tree, Point &w, int k, Point &q);
    static Point_V rtkmethod(RStarTree &treeP, RStarTree &treeW, int k, Point &q);
};


inline double RTK::dot(Point &a, Point &b)
{
    double score = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        score += a.coords[dim] * b.coords[dim];
    return score;
}

inline double RTK::dotMbrLow(Point &point, Mbr &mbr)
{
    double score = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        score += point.coords[dim] * mbr.coord[dim][0];
    return score;
}

inline double RTK::dotMbrUp(Point &point, Mbr &mbr)
{
    double score = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        score += point.coords[dim] * mbr.coord[dim][1];
    return score;
}

inline double RTK::dotMMLow(Mbr &a, Mbr &b)
{
    double score = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        score += a.coord[dim][0] * b.coord[dim][0];
    return score;
}

inline double RTK::dotMMUp(Mbr &a, Mbr &b)
{
    double score = 0;
    for(size_t dim = 0; dim < DIM; dim++)
        score += a.coord[dim][1] * b.coord[dim][1];
    return score;
}





}

#endif /*_DYY_RTK_METHOD_HPP_*/
