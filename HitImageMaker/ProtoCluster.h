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
    : _polygon(points)
  { _hlines.reserve(50); }

  // ProtoCluster(ProtoCluster& p)
  // {
  //   this->_polygon = p._polygon;
  //   //this._hlines  = p._hlines;
  // }
  
  /// Default destructor
  virtual ~ProtoCluster(){}

  Polygon* polygon() { return &_polygon; }

  void AddLine    (float* line)                 { _hlines.push_back(line); }
  // void AddHits(std::vector<larlite::hit>& hits) { _hits = hits;            }

  std::vector<std::pair<float,float> >& hits() { return _hits; }

  //why can't this be reference?
  void add_hit(std::pair<float,float>  h) { return _hits.emplace_back(h); }
  
  // Based on the hough lines that intersect me, determine the most probable direction of myself (with some spread)
  // void DetermineProbableDirection();
  
private:
  
  Polygon _polygon;

  std::vector<float* > _hlines;
  //std::vector<size_t>  _hitidx;
  std::vector<std::pair<float,float> >  _hits;
  

};

#endif


