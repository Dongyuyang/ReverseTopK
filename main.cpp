#include "dyy_RTK_method.hpp"
#include "dyy_data.hpp"
#include <iostream>


typedef dyy::Data DD;

int main()
{
    DD data;

    /*Read data from files*/
    DD::loadPoint("P.data", data.Products);
    DD::loadPoint("W.data", data.Weights);

    /*Build Rtree P and W*/
    DD::buildTree(data.Products, data.entriesP, &data.RtreeP);
    DD::buildTree(data.Weights, data.entriesW, &data.RtreeW);

    dyy::Point q(0.1,0.1);

    /*Reverse top-k query*/
    /*return a vector of w's*/
    // auto rtk = dyy::RTK::rtkmethod(data.RtreeP, data.RtreeW, 10, q );

    /*Reverse k-rank query*/
    /*return w id and rank*/
    auto rkr = dyy::RTK::rkrmethod(data.RtreeP, data.RtreeW,q,10);
    dyy::RTK::print_map(rkr);

    return 1;
}
