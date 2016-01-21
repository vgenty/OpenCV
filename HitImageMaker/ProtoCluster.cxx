#ifndef PROTOCLUSTER_CXX
#define PROTOCLUSTER_CXX

#include "ProtoCluster.h"
#include <numeric>
void ProtoCluster::ComputeDirection() {

  if ( _hlines.size() == 0 )
    return;
  
  _avg_length   = 0;
  _avg_angle    = 0;
  _w_avg_angle  = 0;

  std::vector<float> lengths; lengths.resize(_hlines.size());
  
  for(auto& line : _hlines) {

    auto& x1 = line[0];
    auto& y1 = line[1];
    auto& x2 = line[2];
    auto& y2 = line[3];

    auto dx = x2 - x1;
    auto dy = y2 - y1;

    float pi    = 3.14159;
    float angle = std::atan( dy / dx );

    if ( dy < 0 && dx > 0 ) angle += 2.0 * pi;
    if ( dy > 0 && dx < 0 ) angle += pi;
    if ( dy < 0 && dx < 0 ) angle += pi;
    
    float length = std::sqrt( dx*dx + dy*dy );
    
    _avg_length  += length;
    _avg_angle   += angle;
    _w_avg_angle += angle*length;
    lengths.push_back(length);
  }
  

  _avg_length /= _hlines.size();
  _avg_angle  /= _hlines.size();

  float total_length = std::accumulate(lengths.begin(), lengths.end(), 0);
  _w_avg_angle /= total_length;
  
}




#endif
