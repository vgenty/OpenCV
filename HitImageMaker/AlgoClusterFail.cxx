#ifndef ALGOCLUSTERFAIL_CXX
#define ALGOCLUSTERFAIL_CXX

#include "AlgoClusterFail.h"


namespace larlite {

  AlgoClusterFail::AlgoClusterFail() {
    _dilation_size = 5;
    
    _gauss_blur_size = 5;
    _gauss_sigma_X = 75;

    _thresh = 5; 
    _maxval = 255; 

    _hough_rho = 20.0;

    //_hough_theta;

    _hough_threshold = 100;
    _hough_min_line_length = 100;
    _hough_max_line_gap = 60;  
		 
    _canny_threshold1 = 0;
    _canny_threshold2 = 1;
    _canny_app_size = 3;  


    _merge_min_distance = 20.0;
    _merge_min_angle    = 0.35;// ? right?
    
    import_array();

    init();
  }
  AlgoClusterFail::AlgoClusterFail(const ::fcllite::PSet &pset) {

    _dilation_size   = pset.get<int>("DilationSize");
    
    _gauss_blur_size = pset.get<int>("GausBlurSize");
    _gauss_sigma_X   = pset.get<int>("GausSigmaX");;

    _thresh       = pset.get<int>   ("Thresh");
    _maxval       = pset.get<double>("MaxVal");

    _hough_rho             = pset.get<double>("HoughRho");
    // _hough_theta           = pset.get<double>("HoughTheta");
    _hough_threshold       = pset.get<int>   ("HoughThreshold");
    _hough_min_line_length = pset.get<double>("HoughMinLineLength");
    _hough_max_line_gap    = pset.get<double>("HoughMaxLineGap");
		 
    _canny_threshold1 = pset.get<double>("CannyThreshold1");
    _canny_threshold2 = pset.get<double>("CannyThreshold2");
    _canny_app_size   = pset.get<int>   ("CannyAppSize");

    _merge_min_distance = pset.get<float>("MergeMinDistance");
    _merge_min_angle    = pset.get<float>("MergeMinAngle");
    
    import_array();

    init();

  }
  
  void AlgoClusterFail::init() {
    
    _dilated_v.resize(3);
    _blur_v   .resize(3);
    _binary_v .resize(3);
    _canny_v  .resize(3);

    _houghs_v.resize(3);
    _real_hough_v.resize(3);
    _hulls_v .resize(3);
    
    _possiblebreak_v.resize(3);
    
    _p_clusters_v.resize(3);
    _other_hits_v.resize(3);
  }

  void AlgoClusterFail::reset(const ::cv::Mat& image,size_t plane) {

    auto size = image.size();
    auto type = image.type();

    _dilated_v[plane].release();
    _blur_v   [plane].release();
    _binary_v [plane].release();
    _canny_v  [plane].release();

    std::cout << "in " << __FUNCTION__ << " size : " << size << " type: " << type << "\n";
    _dilated_v[plane].create(size,type);
    _blur_v   [plane].create(size,type);
    _binary_v [plane].create(size,type);
    _canny_v  [plane].create(size,type);

    _houghs_v    [plane].clear();
    _real_hough_v[plane].clear();
    _hulls_v     [plane].clear();

    _possiblebreak_v[plane].clear();
    _p_clusters_v[plane].clear();
    _other_hits_v[plane].clear();
    
  }
  
  void AlgoClusterFail::DecideClusters(event_hit* hits,
				       event_cluster* clusters,
				       AssSet_t* my_ass,
				       const std::vector<::cv::Mat>& images) {
    
    //First thing, we should partition the hits up into planes
    std::map<size_t, std::vector<const hit*> > plane_hits;

    //three planes, reserve the space
    for( unsigned i = 0; i < 3; ++i )

      plane_hits[i].reserve( hits->size() );

      
    //put references to the hits in this map
    for(const auto& hit : *hits)

      plane_hits[hit.WireID().Plane].push_back( &hit );
    
    
    for(const auto& p_hit: plane_hits) { //loop over each plane
      std::cout << "\n\n\n\n\t++++++++++a\n\n\n\n" ;

      auto& p_plane = p_hit.first;
      std::cout << p_plane << "\n";

      auto& p_hits  = p_hit.second;

      if ( p_plane > images.size() - 1 )
	continue;
	
      auto& image   = images.at( p_plane );

      if ( image.size() == cv::Size(0,0) )
	continue;
      
      reset(image,p_plane);
     
      auto& _dilated = _dilated_v[p_plane];
      auto& _blur    = _blur_v   [p_plane];
      auto& _binary  = _binary_v [p_plane];
      // auto& _canny   = _canny_v  [p_plane];

      //auto& _houghs   = _houghs_v    [p_plane];
      auto& _hulls    = _hulls_v     [p_plane];

      auto& _p_clusters = _p_clusters_v[p_plane];

      auto& _other_hits = _other_hits_v[p_plane];

      auto& _possiblebreak = _possiblebreak_v[p_plane];

      //Dilate
      // auto kernel = ::cv::getStructuringElement( ::cv::MORPH_RECT, ::cv::Size(_dilation_size,_dilation_size) );
      // ::cv::dilate(image,_dilated,kernel);
      //Dilate

      auto kernel = ::cv::getStructuringElement( ::cv::MORPH_RECT, ::cv::Size(_dilation_size,_dilation_size) );
      ::cv::dilate(image,_dilated,kernel);

      //Gaussian Blur
      ::cv::GaussianBlur(_dilated,_blur,::cv::Size(_gauss_blur_size,_gauss_blur_size),_gauss_sigma_X);
      
      //Threshold
      //double threshold(InputArray src, OutputArray dst, double thresh, double maxval, int 
      auto what = threshold(_blur,_binary,_thresh,_maxval,::cv::THRESH_BINARY); // what is the return of this?
      
      std::vector<std::vector<cv::Point> > cv_contour_v;
      std::vector<::cv::Vec4i> cv_hierarchy_v;
      ::cv::findContours(_binary,cv_contour_v,cv_hierarchy_v,
    			 CV_RETR_EXTERNAL,
    			 CV_CHAIN_APPROX_SIMPLE);
      
      
      
      //we have the contours, their might be "close" contours where two points are within
      //some minimum distance to each other

      //lets make protoclusters here
      _p_clusters.reserve(cv_contour_v.size());
      
      for(auto& contour : cv_contour_v) {
	
	std::vector<std::pair<float,float> > cont; cont.reserve(contour.size());
	
	for(auto& p : contour)
	  
	  cont.emplace_back(p.x,p.y);
	
	_p_clusters.emplace_back(cont);
	
      }
      

	
      std::map<size_t,std::vector<size_t> > compatable;

      // Are two pclusters compatable
      for(unsigned i = 0; i < _p_clusters.size(); ++i ) {
	
	auto& p1 = _p_clusters[i];

	for(unsigned j = 0; j < _p_clusters.size(); ++j ) {
	    
	  if ( j == i ) continue;
	  
	  auto& p2 = _p_clusters[j];

	  // already connected
	  
	  if ( std::find(compatable[j].begin(), compatable[j].end(), i) != compatable[j].end() )
	    continue;
	  
	  //are they close enough?
	  if (  p1.polygon()->Distance( * p2.polygon() ) > _merge_min_distance  )
	    continue;

	  std::cout << "have distance... " << p1.polygon()->Distance( * p2.polygon() ) << "\n";
	  auto s = p1.polygon()->TwoClosest( *p2.polygon() );

	  std::cout << "s* = " << s.first << "," << s.second << "\n";
	  p1.polygon()->RemoveVertex(s.first);
	  p2.polygon()->RemoveVertex(s.second);
	  
	  compatable[i].push_back(j);

	}
      }
	
      combine_clusters(compatable,_p_clusters);
      convert_convexhull(_p_clusters);

      auto rr = resolve_overlaps(_p_clusters);


      //We need to separate clusters based on output of hough, but how?
      
      //that was terrible.
      std::cout << "\t==> _p_Clusters size...: " << _p_clusters.size() << "\n";

      //write the final hulls to be written to output
      for(unsigned k = 0; k < _p_clusters.size(); ++k) {
	
	_hulls.resize( _p_clusters.size()) ;
	
	for(unsigned p = 0; p < _p_clusters[k].polygon()->Size(); ++p)  {

	  _hulls[k].emplace_back(_p_clusters[k].polygon()->Point(p));

	}
      }
	
      

      //now do the actual clustering, just do conversion locally for now, later we
      //maybe put into utility function
        
      std::vector<int> x_min_v;
      std::vector<int> x_max_v;
      std::vector<int> y_min_v;
      std::vector<int> y_max_v;
      
      std::map<size_t,bool> used_hit;

      for(unsigned i = 0; i < hits->size(); ++i)  {

	auto& h = hits->at(i);
	
	size_t plane = h.WireID().Plane;
	
	if(h.Integral() < 5. || plane != p_plane) {
	  used_hit[i] = true;
	  continue;
	}
	

	if(plane >= x_min_v.size()) {
	  x_min_v.resize(plane+1,4000);
	  x_max_v.resize(plane+1,0);
	  y_min_v.resize(plane+1,9600);
	  y_max_v.resize(plane+1,0);
	}
	
	int wire = h.WireID().Wire;
	int time = (int)(h.PeakTime());
	
	if( x_min_v[plane] > wire ) x_min_v[plane] = wire;
	if( x_max_v[plane] < wire ) x_max_v[plane] = wire;
	if( y_min_v[plane] > time ) y_min_v[plane] = time;
	if( y_max_v[plane] < time ) y_max_v[plane] = time;
	
	used_hit[i] = false;
	
      }

      std::vector<size_t> hit_idx; hit_idx.reserve(p_hits.size());
      
      geo::PlaneID pID(0,0,p_plane);

      int q = 0;
      
      for(auto& c : _p_clusters) {

	hit_idx.clear();
	
	for(size_t i = 0; i < hits->size(); ++i) {

	  if( used_hit[i] ) continue;
	  
	  auto& hit = hits->at(i);

	  float wire     =  (float) hit.WireID().Wire;
	  float time     =  (float) hit.PeakTime();
	  size_t plane   =  hit.WireID().Plane;
	  
	  wire -= (float) x_min_v[plane];
	  time -= (float) y_min_v[plane];
	  
	  if( c.polygon()->PointInside(std::make_pair(time,wire)) )
	    { hit_idx.push_back(i); used_hit[i] = true; c.add_hit( std::make_pair(time,wire) ); }
	  
	  
	}


	if ( hit_idx.size() == 0 )
	  continue;
	
	int ww = 0;

	// hit_idx holds the hit index from all hits. Lets make a new image with hits only in this clusters
	// then run hough on it and see if we can break them up
	
	float xmin = 99999.9;
	float ymin = 99999.9;
	float xmax = 0.0;
	float ymax = 0.0;
	
	
	for( auto& idx : hit_idx ) {

	  auto& hit = hits->at(idx);
	  
	  float wire  =  (float) hit.WireID().Wire;
	  float time  =  (float) hit.PeakTime();

	  if ( wire < xmin ) xmin = wire;
	  if ( time < ymin ) ymin = time;
	  if ( wire > xmax ) xmax = wire;
	  if ( time > ymax ) ymax = time;
	  
	}
	
	::cv::Mat ho(xmax - xmin + 1,
		     ymax - ymin + 1,
		     CV_8UC1, cvScalar(0.));
	  
	for( auto& idx : hit_idx ) {

	  auto& hit = hits->at(idx);

	  float wire     =  (float) hit.WireID().Wire;
	  float time     =  (float) hit.PeakTime();

	  wire -= xmin;
	  time -= ymin;

	  //dont use charge
	  ho.at<unsigned char>(wire,time) = 255;

	}
	
	
	//We can blur it, decide if it's track or shower at this point....
	auto kernel1 = ::cv::getStructuringElement(cv::MORPH_ELLIPSE,::cv::Size(5,5));
	::cv::dilate(ho,ho,kernel1,::cv::Point(-1,-1),2);
	::cv::blur(ho,ho,::cv::Size(5,5));
	auto what2 = ::cv::threshold(ho,ho,0,255,CV_THRESH_BINARY);

	cv_contour_v.clear();
	cv_hierarchy_v.clear();
	::cv::findContours(ho,cv_contour_v,cv_hierarchy_v,
			   CV_RETR_EXTERNAL,
			   CV_CHAIN_APPROX_SIMPLE);

	std::vector<std::vector<::cv::Point> > ctors;
	ctors.resize(cv_contour_v.size());

	for( size_t k = 0; k < cv_contour_v.size(); k++ )
	  ::cv::approxPolyDP(cv_contour_v[k], ctors[k], 0.1, true);
	  
	std::cout << " I found ..." << ctors.size() << " number of contours in ho\n";

	auto largest = int{0};
	size_t index = 0;
	for ( unsigned i = 0; i < ctors.size(); ++i) {
	  std::cout << "\t==>ctors[i].size() = " << ctors[i].size() << "\n";
	  if ( ctors[i].size() > largest ) {
	    largest = ctors[i].size();
	    index = i;
	  }
	}
	
	std::cout << "one with the largest size is: " << index << "\n";

	std::vector<cv::Vec4i> defects;
	std::vector<float> d_d;

	if ( ctors.size() ) {

	  std::vector<int> hullpts;
	  ::cv::convexHull(ctors[index],hullpts);
	  std::cout << "hullpts.size()  " << hullpts.size() << "\n";
	  hullpts.push_back(hullpts.at(0));
	
	  std::cout << "\n(";
	  for ( auto& pt : hullpts ) std::cout << pt << ",";
	  std::cout << ")\n";
	
	  //we use cluster set_width to 0 for shower like 1 for track like ?
	  //
	  
	  ::cv::convexityDefects(ctors[index],hullpts,defects);
	  for(auto & defect : defects) {
	    std::cout << "defect distance = " << defect[3]/256.0 << "\n";
	    d_d.push_back(defect[3]/256.0);
	  }

	  
	}
	
	_possiblebreak.emplace_back(ho);
		
	q++;
	std::cout << "\t--> about to emplace_back\n";
	if( hit_idx.size() > 0 ) {
	  cluster cc;
	  cc.set_id(clusters->size());
	  cc.set_original_producer(clusters->name());
	  cc.set_planeID(pID);
	  cc.set_view(geo::View_t(p_plane));
	  cc.set_defects(d_d);
	  
	  clusters->emplace_back(cc);
	  
	  AssUnit_t one_ass; one_ass.reserve( p_hits.size() );
	  my_ass->emplace_back(one_ass);
	  
	} else { continue; }

	
	auto& this_ass = my_ass->back();

	for(const auto& h : hit_idx)
	  { this_ass.push_back(h); }
	std::cout << "\t--> done emplace_back\n";
	hit_idx.clear();
	std::cout << "\t~~> hit_idx cleared...";

      }

      std::cout << "\t##> this is only temporary\n";
      
      //this is only temporary!!!!
      for(size_t i = 0; i < hits->size(); ++i) {
	  
	if( used_hit[i] ) continue;
	
	auto& hit = hits->at(i);
	
	float wire     =  (float) hit.WireID().Wire;
	float time     =  (float) hit.PeakTime();
	size_t plane   =  hit.WireID().Plane;
	
	wire -= (float) x_min_v[plane];
	time -= (float) y_min_v[plane];

	_other_hits.emplace_back(std::make_pair(time,wire));
	
      }
      
      
    }
    

  }



  int AlgoClusterFail::resolve_overlaps(std::vector<ProtoCluster>& _p_clusters) {

    std::vector<ProtoCluster> combined;
    combined.reserve(_p_clusters.size());

    //Overlaps, holds indicies
    std::vector<size_t> overlaps; overlaps.reserve(_p_clusters.size());
      
    int n = 0; // number of overlaps seen
    int c = 0; // loop counter

    while(true) {

      n = 0;
      std::map<size_t,bool> used;

      if ( c  >  0 ) {
	std::swap(_p_clusters,combined); // should use swap?
	combined.clear();
      }

	
      for(unsigned i = 0; i < _p_clusters.size(); ++i) {
	used[i] = false;
	//is this necessary?
	// _p_clusters[i].polygon()->UntanglePolygon();
      }
      
      for(unsigned c1 = 0; c1 < _p_clusters.size(); ++c1) {

	auto& p_cluster1 = _p_clusters[c1];
	
	overlaps.clear();
	
	if ( used[c1] ) continue;
		  
	for(unsigned c2 = 0; c2 < _p_clusters.size(); ++c2) {
	  
	  if ( c1 == c2 )  continue;
	  if ( used[c2] )  continue;

	  auto& p_cluster2 = _p_clusters[c2];
	  
	  //check if they overlap
	  if ( p_cluster1.polygon()->PolyOverlapSegments( *( p_cluster2.polygon() ) ) ) { 
	    overlaps.push_back(c2);
	    n++;
	  }
	    
	}

	// no overlaps, continue
	if(overlaps.size() == 0)
	  continue;

	//add the first one on
	overlaps.push_back(c1);

	//count the number of points
	int n_points = 0;
	for(const auto& idx : overlaps)

	  n_points += _p_clusters[idx].polygon()->Size();

	//create convex hull out of combined points
	auto out = convert_fewconvex(overlaps,_p_clusters);
	
	//create new protocluster from overlap
	combined.emplace_back(out);
	
	//mark the used overlaps
	for(const auto& idx : overlaps)
	  used[idx] = true;

      } // end c1 loop

	//now add on the clusters that were not combined
      for(const auto& u : used) {
	if ( u.second ) continue;

	// again, can be avoided if we know std::move
	std::vector<std::pair<float,float> > out ( _p_clusters[u.first].polygon()->Size() );

	//loop over points
	for(unsigned k = 0; k < _p_clusters[u.first].polygon()->Size(); ++k) 

	  out[k] = _p_clusters[u.first].polygon()->Point(k);
	
	combined.emplace_back(out);
      }
	
      //did we combine any? if not we are done
      if( n == 0 ) break;

      c++;
	
    }

    std::swap(_p_clusters,combined);
	  
    return c;

  }
    
  void AlgoClusterFail::combine_clusters(std::map<size_t,std::vector<size_t> >& to_combine,
					 std::vector<ProtoCluster>& _p_clusters) {
    

    std::vector<ProtoCluster> combined;
    combined.reserve(_p_clusters.size());

    std::vector<std::pair<float,float> > out; out.reserve(_p_clusters.size());

    std::map<size_t,bool> used;
    for(size_t k = 0; k < _p_clusters.size(); ++k) used[k] = false;
    
    for ( const auto& c_index : to_combine ) {

      if ( c_index.second.size() == 0 ) continue;
      
      std::cout << "proto cluster index: " << c_index.first << " connected to " << c_index.second.size() << " other points! {";
      for(const auto& connected : c_index.second) { std::cout << connected << ","; } std::cout << "}\n";
      out.clear();

      //get the cluster
      auto& clus1 = _p_clusters.at(c_index.first);

      used[ c_index.first ] = true;
      
      //put points inside out
      for(unsigned p = 0; p < clus1.polygon()->Size(); ++p)
	out.emplace_back( clus1.polygon()->Point(p) );

      // loop over the ones that are connected
      for(const auto& connected : c_index.second) {

	//get the cluster
	auto& clus2 = _p_clusters.at(connected);
	used[connected] = true;
	
	//put points inside out
	for(unsigned p = 0; p < clus2.polygon()->Size(); ++p)
	  out.emplace_back( clus2.polygon()->Point(p) );
	
      }
      
      combined.emplace_back(out);
      
    }
    
    
    //now add on the clusters that were not combined
    for(const auto& u : used) {
      if ( u.second ) continue;

      // again, can be avoided if we know std::move
      std::vector<std::pair<float,float> > out2 ( _p_clusters[u.first].polygon()->Size() );

      //loop over points
      for(unsigned k = 0; k < _p_clusters[u.first].polygon()->Size(); ++k) 

	out2[k] = _p_clusters[u.first].polygon()->Point(k);
	
      combined.emplace_back(out2);
    }

    // for(auto& com :combined) { com.polygon()->UntanglePolygon();}

    std::swap(_p_clusters,combined);    
  }


  
  std::vector<std::pair<float,float> > AlgoClusterFail::convert_fewconvex(std::vector<size_t>& idx,
									  std::vector<ProtoCluster>& _pcluster) {


    std::vector<cv::Point> combine; combine.reserve ( 100 );
    std::vector<cv::Point> hul;     hul.reserve     ( 100 );
    
    std::vector<std::pair<float,float> > out; 

    for(const auto& id : idx) {

      for(unsigned p = 0; p < _pcluster[id].polygon()->Size(); ++p) {
	auto point = _pcluster[id].polygon()->Point(p);
	combine.emplace_back(point.first,point.second);
      }
      
    }
    
    convexHull( ::cv::Mat(combine), hul, false );
    out.resize( hul.size() );
    
    for(unsigned k = 0; k < hul.size(); ++k) {
      
      float x = (float) hul[k].x;
      float y = (float) hul[k].y;
      
      out[k] = std::make_pair(x,y);
      
    }
    
    return out;
    
  }




  
  std::vector<std::pair<float,float> > AlgoClusterFail::convert_singleconvex(ProtoCluster& _pcluster) {

    std::vector<cv::Point> combine; combine.reserve ( 100 );
    std::vector<cv::Point> hul;     hul.reserve     ( 100 );
    
    std::vector<std::pair<float,float> > out; 
    
    for(unsigned p = 0; p < _pcluster.polygon()->Size(); ++p) {
      auto point = _pcluster.polygon()->Point(p);
      combine.emplace_back(point.first,point.second);
    }
    
    convexHull( ::cv::Mat(combine), hul, false );
    out.resize( hul.size() );
    
    for(unsigned k = 0; k < hul.size(); ++k) {
      
      float x = (float) hul[k].x;
      float y = (float) hul[k].y;
      
      out[k] = std::make_pair(x,y);
      
    }
    
    return out;
    
  }
  
  void AlgoClusterFail::convert_convexhull(std::vector<ProtoCluster>& _p_clusters) {

    std::vector<ProtoCluster> hulled; hulled.reserve(_p_clusters.size());

    for(unsigned g = 0 ; g < _p_clusters.size() ; ++g)

      hulled.emplace_back(convert_singleconvex(_p_clusters[g]));

    std::swap(_p_clusters,hulled);
  }
  
}
#endif
