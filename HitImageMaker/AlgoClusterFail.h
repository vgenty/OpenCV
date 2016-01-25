//by vic
//vgenty@nevis.columbia.edu

#ifndef ALGOCLUSTERFAIL_H
#define ALGOCLUSTERFAIL_H

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

#include "FhiclLite/PSet.h"

namespace larlite {

  class AlgoClusterFail : public BaseAlgoCluster {

  public:

    /// Default constructor
    AlgoClusterFail();

    ///Alternative ctor
    AlgoClusterFail(const ::fcllite::PSet &pset);
    
    /// Default destructor
    virtual ~AlgoClusterFail(){}

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

    std::vector<std::vector<std::array<float,4> > >& real_hough_v(size_t plane) { return _real_hough_v.at(plane); }

    int possible_break_num(size_t plane) { return _possiblebreak_v[plane].size(); }

    PyObject* possible_break(size_t plane,size_t num) {
      ::larcv::convert::NDArrayConverter converter;
      return converter.toNDArray( _possiblebreak_v[plane][num] );
    }
    
    //temporary
    std::vector<std::pair<float,float> >& other_hits(size_t plane) { return _other_hits_v.at(plane); }
    
  private:
    
    void reset(const ::cv::Mat& image,size_t plane);
    void init ();
      
    std::vector<::cv::Mat> _dilated_v;
    std::vector<::cv::Mat> _blur_v;
    std::vector<::cv::Mat> _binary_v;
    std::vector<::cv::Mat> _canny_v;

    //plane     //possible break
    std::vector<std::vector<::cv::Mat> > _possiblebreak_v;

    std::vector< std::vector<std::vector<std::pair<float,float> > > >_hulls_v;
    std::vector<std::vector<std::array<float,4> > >  _houghs_v;

    std::vector<std::vector<ProtoCluster> > _p_clusters_v;
    
    std::vector< std::vector<std::pair<float,float> > > _other_hits_v;


    int  resolve_overlaps(std::vector<ProtoCluster>& _p_clusters);
    void combine_clusters(std::map<size_t,std::vector<size_t> >& to_combine,
			  std::vector<ProtoCluster>& _p_clusters);

    std::vector<std::pair<float,float> > convert_singleconvex(ProtoCluster& _pcluster);
    std::vector<std::pair<float,float> > convert_fewconvex   (std::vector<size_t>& idx,
							      std::vector<ProtoCluster>& _pcluster);

    //plane     //cluster   //hough
    std::vector<std::vector<std::vector<std::array<float,4> > > > _real_hough_v;
    void convert_convexhull(std::vector<ProtoCluster>& _p_clusters);
    
    
    //Various parameters for image processing set by FHiCL
    int _dilation_size;

    int _gauss_blur_size;
    int _gauss_sigma_X;
    int _gauss_sigma_Y; // unused, sigma_y == sigma_x when unspecified

    int    _thresh;
    double _maxval;

    double _hough_rho;
    double _hough_theta;
    int    _hough_threshold;
    double _hough_min_line_length;
    double _hough_max_line_gap;

    double _canny_threshold1;
    double _canny_threshold2;
    int    _canny_app_size;
    
    float _merge_min_distance;
    float _merge_min_angle;
    
    
  };
}
#endif

