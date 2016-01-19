//by vic
//vgenty@nevis.columbia.edu

#ifndef ALGOIMAGEMAKER_H
#define ALGOIMAGEMAKER_H

#include "BaseAlgoImage.h"

struct _object;
typedef _object PyObject;

#ifndef __CINT__
#include "Python.h"
#include "numpy/arrayobject.h"
#endif

namespace larlite {

  class AlgoImageMaker : public BaseAlgoImage {

  public:

    /// Default constructor
    AlgoImageMaker() : _name_s("AlgoImageMaker") {}

    /// Default destructor
    virtual ~AlgoImageMaker(){}

    void CreateImage(const event_hit* ev_hit);
    const ::cv::Mat* GetImage(const size_t plane) { return &_mat_v.at(plane); }    
    
    PyObject* GetPyImage(const size_t plane);


    
  private:

    void clean();
    
    std::string _name_s;
  
    std::vector<int>   _x_min_v;
    std::vector<int>   _x_max_v;
    std::vector<int>   _y_min_v;
    std::vector<int>   _y_max_v;
    std::vector<float> _q_max_v;
  
  };
}

#endif


