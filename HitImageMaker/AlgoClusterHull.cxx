#ifndef ALGOCLUSTERHULL_CXX
#define ALGOCLUSTERHULL_CXX

#include "AlgoClusterHull.h"



namespace larlite {

  AlgoClusterHull::AlgoClusterHull() {
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
    
    import_array();

    init();
  }
  AlgoClusterHull::AlgoClusterHull(const ::fcllite::PSet &pset) {

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

    import_array();

    init();

  }
  
  void AlgoClusterHull::init() {
    
    _dilated_v.resize(3);
    _blur_v   .resize(3);
    _binary_v .resize(3);
    _canny_v  .resize(3);

    _houghs_v.resize(3);
    _hulls_v .resize(3);
    _contours_v.resize(3);
    _p_clusters_v.resize(3);
    _other_hits_v.resize(3);
  }

  void AlgoClusterHull::reset(const ::cv::Mat& image,size_t plane) {

    auto size = image.size();
    auto type = image.type();

    _dilated_v[plane].release();
    _blur_v   [plane].release();
    _binary_v [plane].release();
    _canny_v  [plane].release();
    
    
    _dilated_v[plane].create(size,type);
    _blur_v   [plane].create(size,type);
    _binary_v [plane].create(size,type);
    _canny_v  [plane].create(size,type);

    _contours_v[plane].clear();

    _houghs_v [plane].clear();
    _hulls_v  [plane].clear();

    _p_clusters_v[plane].clear();
    _other_hits_v[plane].clear();
    
  }
  
  void AlgoClusterHull::DecideClusters(event_hit* hits,
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

      auto& p_plane = p_hit.first;
      auto& p_hits  = p_hit.second;
      auto& image   = images.at( p_plane );

      reset(image,p_plane);
     
      auto& _dilated = _dilated_v[p_plane];
      auto& _blur    = _blur_v   [p_plane];
      auto& _binary  = _binary_v [p_plane];
      auto& _canny   = _canny_v  [p_plane];

      auto& _houghs   = _houghs_v[p_plane];
      auto& _hulls    = _hulls_v [p_plane];

      auto& _p_clusters = _p_clusters_v[p_plane];
      auto& _other_hits = _other_hits_v[p_plane];
      auto& _contours   = _contours_v[p_plane];
	
      //Dilate
      auto kernel = ::cv::getStructuringElement( ::cv::MORPH_RECT, ::cv::Size(_dilation_size,_dilation_size) );
      ::cv::dilate(image,_dilated,kernel);

      //Gaussian Blur
      ::cv::GaussianBlur(_dilated,_blur,::cv::Size(_gauss_blur_size,_gauss_blur_size),_gauss_sigma_X);
      
      //Threshold
      //double threshold(InputArray src, OutputArray dst, double thresh, double maxval, int 
      auto what = threshold(_blur,_binary,_thresh,_maxval,::cv::THRESH_BINARY); // what is the return of this?
      
      //HoughLinesP
      //void HoughLinesP(InputArray image, OutputArray lines,
      //                 double rho, double theta, int threshold,
      //                 double minLineLength=0, double maxLineGap=0 )
      
      auto pi = double{3.14159};
      std::vector<::cv::Vec4i> lines;
      ::cv::HoughLinesP(_binary, lines, _hough_rho, pi/180.0,
			_hough_threshold, _hough_min_line_length, _hough_max_line_gap);
      
      _houghs.resize(lines.size());
      for(unsigned i = 0; i < lines.size(); ++i )
	_houghs[i] = { (float) lines[i][0], (float) lines[i][1],
		       (float) lines[i][2], (float) lines[i][3] };      
      
      //Canny
      ::cv::Canny(_binary,_canny,_canny_threshold1,_canny_threshold2,_canny_app_size);
      
      // ==> page 100 of ``Computer Vision with OpenCV"
      // In OpenCV, each individual contour is stored as a vector of points,
      // and all the contours are stored as a vector of contours (i.e. a vector of vectors of points).

      //Contours
      std::vector<std::vector<cv::Point> > cv_contour_v;
      std::vector<::cv::Vec4i> cv_hierarchy_v;
      ::cv::findContours(_binary,cv_contour_v,cv_hierarchy_v,
      			 CV_RETR_EXTERNAL,
      			 CV_CHAIN_APPROX_SIMPLE);
      
      // ::cv::findContours(_canny,cv_contour_v,cv_hierarchy_v,
      // 			 CV_RETR_EXTERNAL,
      // 			 CV_CHAIN_APPROX_NONE);
      
      /////////////////////////////////////////////////////
      //This may be logical break here for new framework?
      //
      //We can start building the skeleton clusters....
      //
      
      
      // convexHull first, fill the protoclusters with convex hull'd contours

      _contours.reserve( cv_contour_v.size() );

      for(auto& contour : cv_contour_v) {
	std::vector<std::pair<float,float> > _contour; _contour.reserve( contour.size() );
	for(auto& point : contour) {
	  _contour.push_back(std::make_pair(point.x,point.y));
	}
	_contours.emplace_back(_contour);
      }
	  
      
      std::vector<std::vector<::cv::Point> > hull( cv_contour_v.size() );

      std::vector<std::vector<std::pair<float,float> > > hulls;
      hulls.resize(cv_contour_v.size());


      _p_clusters.reserve( cv_contour_v.size() );

      for( unsigned i = 0; i < cv_contour_v.size(); i++ ) {

	convexHull( ::cv::Mat(cv_contour_v[i]), hull[i], false );

	hulls[i].resize(hull[i].size());

	for(unsigned k = 0; k < hull[i].size(); ++k) {
	  float x = (float) hull[i][k].x;
	  float y = (float) hull[i][k].y;

	  hulls[i][k] = std::make_pair(x,y);
	}
	
	_p_clusters.emplace_back(hulls[i]);
      }


      //For now just combine overlapping clusters (polygons) can we do recursively?
      std::vector<ProtoCluster> combined;
      combined.reserve(cv_contour_v.size());

      //Overlaps, holds indicies
      std::vector<size_t> overlaps; overlaps.reserve(_p_clusters.size());
      
      int n = 0; // number of overlaps seen
      int c = 0; // loop counter

      while(true) {
	
	n = 0;
	std::map<size_t,bool> used;

	if ( c  >  0 ) {
	  _p_clusters = combined; // should use swap!
	  combined.clear();
	}

	
	for(unsigned i = 0; i < _p_clusters.size(); ++i) {
	  used[i] = false;
	  //is this necessary?
	  _p_clusters[i].polygon()->UntanglePolygon();
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
	  std::vector<cv::Point> combine; combine.reserve ( n_points );
	  std::vector<cv::Point> hul;     hul.reserve     ( n_points );
	
	  std::vector<std::pair<float,float> > out;

	  for(const auto& idx : overlaps) {
	    for(unsigned p = 0; p < _p_clusters[idx].polygon()->Size(); ++p) {
	      auto point = _p_clusters[idx].polygon()->Point(p);
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

      //swap the data
      std::swap(_p_clusters,combined);

      //that was terrible.
      std::cout << "\t==> _P_Clusters size...: " << _p_clusters.size() << "\n";

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
	
	q++;
	// std::cout << "We found ..." << hit_idx.size() << " inside this cluster...\n";
	
	if( hit_idx.size() > 0 ) {
	  cluster cc;
	  cc.set_id(clusters->size());
	  cc.set_original_producer(clusters->name());
	  cc.set_planeID(pID);
	  cc.set_view(geo::View_t(p_plane));
	  
	  clusters->emplace_back(cc);
	  
	  AssUnit_t one_ass; one_ass.reserve( p_hits.size() );
	  my_ass->emplace_back(one_ass);
	  
	} else { continue; }
	
	auto& this_ass = my_ass->back();

	for(const auto& h : hit_idx)
	  { this_ass.push_back(h); }
		      
	hit_idx.clear();
	
      }

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


    

}
#endif
