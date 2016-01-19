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
    
    
    //just do one plane for now
    _algo_image->CreateImage(ev_hit);

    //BaseAlgoImage will draw the pictures for us ( 3 planes );
    //    for( unsigned i = 0; i < 3; ++i ) {
    //auto p = _algo_image->GetImage(2); // get the pointer to plane 2 image
    //}
    
    
    //run the clustering algo (for now give it all the hits and give it all the images
    //later we can break this up
    auto clusters = _algo_cluster->DecideClusters(ev_hit,
						  _algo_image->_mat_v);
    
    return true;
  }
  
  bool HitImageManager::finalize() {
    
    return true;
  }
  

}
#endif
