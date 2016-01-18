#ifndef LARLITE_HITIMAGEMANAGER_CXX
#define LARLITE_HITIMAGEMANAGER_CXX

#include "HitImageManager.h"

namespace larlite {

  bool HitImageManager::initialize() {

    std::cout << "\t==> Initilizing HitImageManager\n";
    
    return true;
  }
  
  bool HitImageManager::analyze(storage_manager* storage) {

    //Manage takes the hits
    auto const ev_hit = storage->get_data<event_hit>("gaushit");  
    
    //BaseAlgoImage will draw the pictures for us ( 3 planes );
    //    for( unsigned i = 0; i < 3; ++i ) {

    //}

    //just do one plane for now
    _algo_image->CreateImage(ev_hit);
    auto p = _algo_image->GetImage(2); // get the pointer to plane 2 image
    
    //for each plane, run the clustering algo
    //auto clusters = _algo_cluster->DecideClusters(ev_hit,p);
    
    return true;
  }
  
  bool HitImageManager::finalize() {
    
    return true;
  }
  

}
#endif
