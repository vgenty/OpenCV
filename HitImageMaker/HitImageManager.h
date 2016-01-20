/**
 * \file HitImageManager.h
 *
 * \ingroup HitImageMaker
 * 
 * \brief Class def header for a class HitImageManager
 *
 * @author vgenty
 */

/** \addtogroup HitImageMaker

    @{*/

#ifndef LARLITE_HITIMAGEMANAGER_H
#define LARLITE_HITIMAGEMANAGER_H

#include "Analysis/ana_base.h"

#include "BaseAlgoImage.h"
#include "BaseAlgoCluster.h"

namespace larlite {

  class HitImageManager : public ana_base{
  
  public:

    /// Default constructor
    HitImageManager(){ _name="HitImageManager"; _fout=0;}

    /// Default destructor
    virtual ~HitImageManager(){}

    virtual bool initialize();
    virtual bool analyze(storage_manager* storage);
    virtual bool finalize();

    void SetImageMaker  (BaseAlgoImage*   algo_image)   { _algo_image   = algo_image;   }
    void SetClusterMaker(BaseAlgoCluster* algo_cluster) { _algo_cluster = algo_cluster; }

    
  protected:
    // ?
    
  private:
    
    BaseAlgoImage*   _algo_image;
    BaseAlgoCluster* _algo_cluster;

  };
}
#endif
