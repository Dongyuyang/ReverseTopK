#include "dyy_data.hpp"




namespace dyy{

/*********************************************************************
 * Point
 ********************************************************************/
std::istream& operator>>(std::istream &in, Point &p)
{
    p.coords.resize(DIM);
    for(size_t d = 0; d < DIM; d++)
        in >> p.coords[d];
    return in;
}

bool Point::operator<(const dyy::Point &p) const
{
    for(size_t dim = 0; dim < DIM; dim++)
        if(coords[dim] != p.coords[dim])
            return coords[dim] < p.coords[dim];
    return true;
}

void Point::print()
{
    std::cout << "<" << coords[0];
    for(size_t dim = 1; dim < DIM; dim++)
        std::cout << "," << coords[dim];
    std::cout << ">" << std::endl;
}

/*********************************************************************
 * Data
 ********************************************************************/

void Data::loadPoint(std::string data_file,  Point_V &points)
{
    points.clear();
    std::ifstream in(data_file.c_str());
    assert(in.is_open());
    int id = 0;
    while(true){
        Point point;
        in >> point;
        if(in){
            point.id = id++;
            points.push_back(point);
        }
        else
            break;
    }
    in.close();
}

void Data::buildTree(Point_V& points, Entry_V& entries, RStarTree* tree)
{
    entries.clear();
    for(size_t ip = 0; ip < points.size(); ip++){
        Data_P datap = &points.at(ip);
        Mbr mbr(points.at(ip).coords);
        LeafNodeEntry entry(mbr, datap);
        entries.push_back(entry);
    }

    std::cout << entries.size() << " entries created" << std::endl;

    for(size_t ie = 0; ie < entries.size(); ie++)
        tree->insertData(&entries.at(ie));

    std::cout << tree->root->aggregate << " entries created" << std::endl;
}



}
