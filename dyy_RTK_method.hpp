#ifndef _DYY_RTK_METHOD_HPP_
#define _DYY_RTK_METHOD_HPP_

#include "dyy_data.hpp"
#include <map>
#include <vector>
#include <iterator>

namespace dyy{


class RTK
{
public:

    /*dot procut*/
    static double dot(Point &a, Point &b);
    static double dotMbrLow(Point &point, Mbr &mbr);
    static double dotMbrUp(Point &point, Mbr &mbr);
    static double dotMMLow(Mbr &a, Mbr &b);
    static double dotMMUp(Mbr &a, Mbr &b);

    /*Reverse top-k method*/
    static bool inkPW(RStarTree &tree, Mbr &Ew, int k, Point &q);
    static bool ink(RStarTree &tree, Point &w, int k, Point &q);
    static Point_V rtkmethod(RStarTree &treeP, RStarTree &treeW, int k, Point &q);

    /*Reverse k-Rank method*/
    class ARankResult
    {
    public:
        bool isBetter;
        int rank;
        /*
          1:All in
          flag =    0:Need check children
          -1:All out
        */
        int flag;
    };
    static ARankResult inRank(RStarTree &tree, Point &w, int minRank, Point &q);
    static ARankResult inRankPW(RStarTree &tree, Mbr &Ew, int minRank, Point &q);
    typedef std::multimap<int, int> BUFFER;
    static void update_map(BUFFER &map, int key, int value);
    static void print_map(BUFFER &map);
    static BUFFER rkrmethod(RStarTree &treeP, RStarTree &treeW, Point &q, int k);

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

inline void RTK::print_map(BUFFER &map)
{
    for ( std::multimap<int,int>::const_iterator it = map.begin();
          it != map.end(); it++)
        {
            std::cout << "ARank: " << it->first;
            std::cout << " ,Id: " << it->second <<std::endl;
        }
}

inline void RTK::update_map(BUFFER &map, int key, int value)
{
    map.insert(std::pair<int,int>(key,value));
    map.erase(std::prev( map.end()));
}





}

#endif /*_DYY_RTK_METHOD_HPP_*/
