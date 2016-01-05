#ifndef CONTOURTREEMAKER_CXX
#define CONTOURTREEMAKER_CXX

#include "ContourTreeMaker.h"

namespace imutil {

  ContourTreeMaker::ContourTreeMaker(const std::string name)
    : _tree(name.c_str(),"Contour TTree")
  {
    _tree.Branch("contours","std::vector<std::vector<std::pair<double,double> > >",&_contour_array);
  }

  void ContourTreeMaker::FillContour(const ContourID_t index,
				     const double x,
				     const double y)
  {
    if(index >= _contour_array.size()) _contour_array.resize(index+1);
    _contour_array[index].emplace_back(x,y);
  }

}

#endif
