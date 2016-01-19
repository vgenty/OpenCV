//by vic
//vgenty@nevis.columbia.edu

#ifndef BASEALGOCLUSTER_H
#define BASEALGOCLUSTER_H

#include <iostream>

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include "opencv2/core/mat.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace larlite {

  class BaseAlgoCluster{

  public:

    /// Default constructor
    BaseAlgoCluster(){}

    /// Default destructor
    virtual ~BaseAlgoCluster(){}

    virtual event_cluster DecideClusters(event_hit* hits,
					 const std::vector<::cv::Mat>& images) = 0;
    
  };

}

#endif

