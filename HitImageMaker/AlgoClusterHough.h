//by vic
//vgenty@nevis.columbia.edu

#ifndef ALGOCLUSTERHOUGH_H
#define ALGOCLUSTERHOUGH_H

#include <vector>
#include <map>
#include <array>

#include "BaseAlgoCluster.h"

#include "Utils/DataTypes.h"

struct _object;
typedef _object PyObject;

#ifndef __CINT__
#include "Python.h"
#include "numpy/arrayobject.h"
#endif

#include "Utils/NDArrayConverter.h"

namespace larlite {

  class AlgoClusterHough : public BaseAlgoCluster {

  public:

    /// Default constructor
    AlgoClusterHough(){ import_array(); }

    /// Default destructor
    virtual ~AlgoClusterHough(){}

    void  DecideClusters(event_hit* hits,
			 event_cluster* clusters,
			 AssSet_t* my_ass,
			 const std::vector<::cv::Mat>& images);
    
 
    PyObject* GetDilatedImage() { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray(_dilated); }
    PyObject* GetBlurImage()    { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray(_blur);    }
    PyObject* GetBinaryImage()  { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray(_binary);  }
    PyObject* GetCannyImage()   { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray(_canny);   }


    
    PyObject* GetContour(const size_t contour_index);
    size_t NumContours() const;
    
    size_t contoursize() { return _contour_v2.size(); }
    std::vector<std::pair<double,double> >& contour(const size_t c) { return _contour_v2.at(c); }

    std::vector< std::vector<std::pair<float,float> > >& hulls()     { return _hulls; }
    std::vector< std::vector<std::pair<float,float> > >& hulls2()    { return _hulls2; }
    std::vector<std::pair<float,float> >               & plane2pts() { return _plane2pts; }


    std::vector<std::array<float,4> >& houghs() { return _houghs; }
    
  private:
    
    void init(const ::cv::Mat& image);
      
    ::cv::Mat _dilated;
    ::cv::Mat _blur;
    ::cv::Mat _binary;
    ::cv::Mat _canny;

    std::vector<larcv::Point2DArray> _contour_v;

    std::vector<std::vector<std::pair<double,double> > > _contour_v2;

    std::vector<std::array<float,4> > _houghs;

    std::vector<std::vector<std::pair<float,float> > > _hulls;

    std::vector<std::vector<std::pair<float,float> > > _hulls2;

    std::vector<std::pair<float,float> > _plane2pts;
    
  };
}
#endif

