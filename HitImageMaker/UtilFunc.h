#ifndef UTILFUNC_H
#define UTILFUNC_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>


#include "TH1D.h"

//#include "OpticalRecoTypes.h"

namespace larlite {
  
  double mean(const std::vector<short>& wf, size_t start=0, size_t nsample=0);

  double edge_aware_mean(const std::vector<short>& wf, int start, int end);
  
  double std(const std::vector<short>& wf, const double ped_mean, size_t start=0, size_t nsample=0);

  double BinnedMaxOccurrence(const std::vector<float>& mean_v,const size_t nbins);

  double BinnedMaxTH1D(const std::vector<float>& v ,int bins);
  
  int sign(double val);  

}

#endif
