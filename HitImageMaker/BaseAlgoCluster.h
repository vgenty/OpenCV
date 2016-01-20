//by vic
//vgenty@nevis.columbia.edu

#ifndef BASEALGOCLUSTER_H
#define BASEALGOCLUSTER_H

#include <iostream>

#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/event_ass.h"

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

    virtual void DecideClusters(event_hit* hits,
				event_cluster* clusters,
				AssSet_t* my_ass,
				const std::vector<::cv::Mat>& images) = 0;
    
  };

}

#endif

