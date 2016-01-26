#ifndef ALGOIMAGEMAKER_CXX
#define ALGOIMAGEMAKER_CXX

#include "AlgoImageMaker.h"
#include "Utils/NDArrayConverter.h"

namespace larlite {

  void AlgoImageMaker::CreateImage(const event_hit* ev_hit) {
    
    clean();
    
    for(auto const& h : *ev_hit) {
      
      if(h.Integral()<5.) continue;

      size_t plane = h.WireID().Plane;
      if(plane >= _x_min_v.size()) {
	_x_min_v.resize(plane+1,4000);
	_x_max_v.resize(plane+1,0);
	_y_min_v.resize(plane+1,9600);
	_y_max_v.resize(plane+1,0);
	_q_max_v.resize(plane+1,0);
      }

      int wire = h.WireID().Wire;
      int time = (int)(h.PeakTime());
      float  q = h.Integral();

      if( _x_min_v[plane] > wire ) _x_min_v[plane] = wire;
      if( _x_max_v[plane] < wire ) _x_max_v[plane] = wire;
      if( _y_min_v[plane] > time ) _y_min_v[plane] = time;
      if( _y_max_v[plane] < time ) _y_max_v[plane] = time;
      if( _q_max_v[plane] < q    ) _q_max_v[plane] = q;
    }

    // std::cout << "x min v size : " << _x_min_v.size() << "\n";
    
    for(size_t plane=0; plane<_x_min_v.size(); ++plane) {
      // std::cout << "plane " << plane << "\n";
      // std::cout << "max,minX" << _x_max_v[plane] - _x_min_v[plane] + 1 << "\n";
      // std::cout << "max,minY" << _y_max_v[plane] - _y_min_v[plane] + 1 << "\n";

      if ( _x_max_v[plane] - _x_min_v[plane] + 1 > 0 && _y_max_v[plane] - _y_min_v[plane] + 1 > 0) {
	::cv::Mat mat(_x_max_v[plane] - _x_min_v[plane] + 1,
		      _y_max_v[plane] - _y_min_v[plane] + 1,
		      CV_8UC1, cvScalar(0.));
	
	_mat_v.emplace_back(mat);
      }
      else
	_mat_v.emplace_back(0,0,CV_8UC1,cvScalar(0.));
      
    }

    // std::cout << "(" << _x_min_v[2] << ","
    // 	      << _y_min_v[2] << ")\n";
      
    
    for(auto const& h : *ev_hit) {

      if(h.Integral()<5.) continue;

      int wire     = h.WireID().Wire;
      int time     = (int)(h.PeakTime());
      size_t plane = h.WireID().Plane;
      int charge   = ((256. * h.Integral() / _q_max_v[plane]));

      wire -= _x_min_v[plane];
      time -= _y_min_v[plane];

      auto& mat = _mat_v[plane];

      charge += mat.at<unsigned char>(wire,time);

      if(charge>256) charge = 256;
      mat.at<unsigned char>(wire,time) = (unsigned char)(charge);
	
    }
    
  }
  
  void AlgoImageMaker::clean() {

    _x_min_v.clear();
    _x_max_v.clear();
    _y_min_v.clear();
    _y_max_v.clear();
    _q_max_v.clear();
    
    _mat_v.clear();
  }


  PyObject* AlgoImageMaker::GetPyImage(const size_t plane)
  {

    if(plane >= _mat_v.size()) {
      // std::cout << "\t==X Plane doesn't exist\n";
      throw std::exception();
    }
    
    ::larcv::convert::NDArrayConverter converter;
    return converter.toNDArray(_mat_v[plane]);

  }
  
}
#endif
