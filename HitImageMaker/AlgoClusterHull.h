//by vic
//vgenty@nevis.columbia.edu

#ifndef ALGOCLUSTERHULL_H
#define ALGOCLUSTERHULL_H

#include <vector>
#include <map>
#include <array>

#include "BaseAlgoCluster.h"

#include "Utils/DataTypes.h"

#include "ProtoCluster.h"

struct _object;
typedef _object PyObject;

#ifndef __CINT__
#include "Python.h"
#include "numpy/arrayobject.h"
#endif

#include "Utils/NDArrayConverter.h"

namespace larlite {

  class AlgoClusterHull : public BaseAlgoCluster {

  public:

    /// Default constructor
    AlgoClusterHull(){ import_array(); init(); }

    /// Default destructor
    virtual ~AlgoClusterHull(){}

    void DecideClusters(event_hit* hits,
			event_cluster* clusters,
			AssSet_t* my_ass,
			const std::vector<::cv::Mat>& images);
    
 
    PyObject* GetDilatedImage(size_t plane) { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray( _dilated_v[plane] ); }
    PyObject* GetBlurImage   (size_t plane) { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray( _blur_v   [plane] ); }
    PyObject* GetBinaryImage (size_t plane) { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray( _binary_v [plane] ); }
    PyObject* GetCannyImage  (size_t plane) { ::larcv::convert::NDArrayConverter converter; return converter.toNDArray( _canny_v  [plane] ); }


    std::vector< std::vector<std::pair<float,float> > >& hulls(size_t plane)  { return _hulls_v.at(plane);   }
    std::vector<std::array<float,4> >& houghs(size_t plane)                   { return _houghs_v.at(plane);  }
    std::vector<ProtoCluster>& p_clusters(size_t plane)                       { return _p_clusters_v.at(plane); }

    //temporary
    std::vector<std::pair<float,float> >& other_hits(size_t plane) { return _other_hits_v.at(plane); }
    
  private:
    
    void reset(const ::cv::Mat& image,size_t plane);
    void init ();
      
    std::vector<::cv::Mat> _dilated_v;
    std::vector<::cv::Mat> _blur_v;
    std::vector<::cv::Mat> _binary_v;
    std::vector<::cv::Mat> _canny_v;


    std::vector< std::vector<std::vector<std::pair<float,float> > > >_hulls_v;
    std::vector<std::vector<std::array<float,4> > >  _houghs_v;

    std::vector<std::vector<ProtoCluster> > _p_clusters_v;
    
    std::vector< std::vector<std::pair<float,float> > > _other_hits_v;
    
  };
}
#endif

