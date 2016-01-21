//by vic
//vgenty@nevis.columbia.edu

#ifndef PROTOCLUSTER_H
#define PROTOCLUSTER_H

#include "Polygon.h"
#include "DataFormat/hit.h"

class ProtoCluster {
  
public:

  /// Default constructor
  ProtoCluster() {}
  
  ProtoCluster(const std::vector< std::pair<float,float> >& points)
    :
    _polygon(points),
    _avg_length(-1),
    _avg_angle(-1),
    _w_avg_angle(-1)
    
  { _hlines.reserve(50); }

  // ProtoCluster(ProtoCluster& p)
  // {
  //   this->_polygon = p._polygon;
  //   //this._hlines  = p._hlines;
  // }
  
  /// Default destructor
  virtual ~ProtoCluster(){}

  Polygon* polygon() { return &_polygon; }

  void AddLine(float* line){ _hlines.push_back(line); }
  void ComputeDirection();
  
  std::vector<std::pair<float,float> >& hits() { return _hits; }

  //why can't this be reference?
  void add_hit(std::pair<float,float>  h) { return _hits.emplace_back(h); }

  float _avg_length;
  float _avg_angle;
  float _w_avg_angle;

  size_t n_lines() { return _hlines.size(); }
    
private:
  
  Polygon _polygon;

  std::vector<float* > _hlines;
  std::vector<std::pair<float,float> >  _hits;
  

};

#endif


