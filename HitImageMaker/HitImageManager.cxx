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
    // std::cout << "creating image.........................\n";
    _algo_image->CreateImage(ev_hit);
    

    //run the clustering algo (for now give it all the hits and give it all the images
    //later we can break this up

    auto ev_cluster = storage->get_data<event_cluster>("hullcluster");
    auto ev_ass     = storage->get_data<event_ass>(ev_cluster->name());

    AssSet_t my_ass;
    // std::cout << "deciding clusters.........................\n";
    _algo_cluster->DecideClusters(ev_hit,
				  ev_cluster,
				  &my_ass,
				  _algo_image->_mat_v);
    
    // std::cout << "setting association.........................\n";
    ev_ass->set_association( ev_cluster->id(), ev_hit->id(), my_ass );
    
    storage->set_id(storage->run_id(),
		    storage->subrun_id(),
		    storage->event_id());
    
    // std::cout << "Found..." << ev_cluster->size() << " clusters \n";
    
    return true;
  }
  
  bool HitImageManager::finalize() {
    
    return true;
  }
  

}
#endif
