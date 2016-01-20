//by vic
//vgenty@nevis.columbia.edu

#ifndef BASEALGOIMAGE_H
#define BASEALGOIMAGE_H

#include <iostream>
#include <string>
#include <vector>

#include "opencv2/core/mat.hpp"

#include "DataFormat/hit.h"

namespace larlite {

  class BaseAlgoImage {

  public:

    /// Default constructor
    BaseAlgoImage() : _name_s("BaseAlgoImage") {}

    /// Default destructor
    virtual ~BaseAlgoImage(){}

    virtual       void       CreateImage(const event_hit* ev_hit) = 0;
    virtual const ::cv::Mat* GetImage   (const size_t plane) = 0;
    
    // PyObject* GetImage(const size_t plane);
    std::vector<::cv::Mat>   _mat_v;
    
  private:
    std::string _name_s;
  
  };
}

#endif
/** @} */ // end of doxygen group 

