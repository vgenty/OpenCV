#ifndef ALGOCLUSTERHOUGH_CXX
#define ALGOCLUSTERHOUGH_CXX

#include "AlgoClusterHough.h"

#include "ProtoCluster.h"

namespace larlite {

  void AlgoClusterHough::init(const ::cv::Mat& image) {

    auto size = image.size();
    auto type = image.type();


    _dilated.release();
    _blur   .release();
    _binary .release();
    _canny  .release();
    
    
    _dilated.create(size,type);
    _blur   .create(size,type);
    _binary .create(size,type);
    _canny  .create(size,type);

    _contour_v.clear();
    _contour_v2.clear();
    _houghs.clear();
    _hulls.clear();
    _hulls2.clear();
  }
  
  event_cluster AlgoClusterHough::DecideClusters(event_hit* hits,
						 const std::vector<::cv::Mat>& images) {
    
    
    event_cluster evt_clusters;
    
    
    //First thing, we should partition the hits up into planes
    std::map<size_t, std::vector<const hit*> > plane_hits;

    //three planes, reserve the space
    for( unsigned i = 0; i < 3; ++i )

      plane_hits[i].reserve( hits->size() );
    
    //put references to the hits in this map
    for(const auto& hit : *hits)
      
      plane_hits[hit.WireID().Plane].push_back( &hit );
    
    
    for(const auto& p_hits: plane_hits) {


      auto& plane = p_hits.first;
      auto& hits  = p_hits.second;
      auto& image = images.at( plane );

      //temp
      if (plane != 2) continue;
      //temp

      
      //Set it u
      init(image);
      
      //Dilate
      std::cout << "\t(" << __FUNCTION__ << ")==> Dilating\n";
      auto kernel = ::cv::getStructuringElement( ::cv::MORPH_RECT, ::cv::Size(5,5) );

      ::cv::dilate(image,_dilated,kernel);

      //Gaussian Blur
      std::cout << "\t(" << __FUNCTION__ << ")==> Bluring\n";
      ::cv::GaussianBlur(_dilated,_blur,::cv::Size(5,5),75);
      
      //Threshold
      //double threshold(InputArray src, OutputArray dst, double thresh, double maxval, int 
      std::cout << "\t(" << __FUNCTION__ << ")==> Threshing\n";
      auto what = threshold(_blur,_binary,5,255,::cv::THRESH_BINARY); // what is the return of this?
      

      //HoughLinesP
      //void HoughLinesP(InputArray image, OutputArray lines,
      //                 double rho, double theta, int threshold,
      //                 double minLineLength=0, double maxLineGap=0 )
      
      auto pi = double{3.14159};
      std::vector<::cv::Vec4i> lines;
      std::cout << "\t(" << __FUNCTION__ << ")==> My friend Hough\n";
      ::cv::HoughLinesP(_binary, lines, 20.0, pi/180.0, 100, 100, 60);
      
      _houghs.resize(lines.size());
      for(unsigned i = 0; i < lines.size(); ++i )
	_houghs[i] = { (float) lines[i][0], (float) lines[i][1],
		       (float) lines[i][2], (float) lines[i][3] };      
      
      //Canny
      std::cout << "\t(" << __FUNCTION__ << ")==> Canny!\n";
      ::cv::Canny(_binary,_canny,0,1,3);
      
      // ==> page 100 of ``Computer Vision with OpenCV"
      // In OpenCV, each individual contour is stored as a vector of points,
      // and all the contours are stored as a vector of contours (i.e. a vector of vectors of points).

      //Contours
      std::cout << "\t(" << __FUNCTION__ << ")==> Contour\n";
      std::vector<std::vector<cv::Point> > cv_contour_v;
      std::vector<::cv::Vec4i> cv_hierarchy_v;
      ::cv::findContours(_canny,cv_contour_v,cv_hierarchy_v,
    			 CV_RETR_EXTERNAL,
    			 CV_CHAIN_APPROX_SIMPLE);
      
      std::cout<< "\t(" << __FUNCTION__ << ")==> Found " << cv_contour_v.size()<<" contours..."<<std::endl;

      
      //
      //We can start building the skeleton clusters.
      //
      
      
      // convexHull first, fill the protoclusters
      std::cout<< "\t(" << __FUNCTION__ << ")==> Building protoclusters\n";
      std::vector<std::vector<::cv::Point> > hull( cv_contour_v.size() );

      //FUCK ME
      _hulls.resize(cv_contour_v.size());

      std::vector<ProtoCluster> p_clusters_v;
      p_clusters_v.reserve( cv_contour_v.size() );

      for( unsigned i = 0; i < cv_contour_v.size(); i++ ) {

	convexHull( ::cv::Mat(cv_contour_v[i]), hull[i], false );

	_hulls[i].resize(hull[i].size());

	for(unsigned k = 0; k < hull[i].size(); ++k) {

	  float x = (float) hull[i][k].x;
	  float y = (float) hull[i][k].y;

	  _hulls[i][k] = std::make_pair(x,y);
	  
					
	}
	
	p_clusters_v.emplace_back(_hulls[i]);

      }


      // TODO : USE HOUGHLINE TO COMBINE PROTOCLUSTERS...

      //for now just combine overlapping clusters (polygons)
      std::vector<ProtoCluster> combined;
      bool more_overlaps = true;
      int n = 0;
      int c = 0;

      while(more_overlaps) {
	
	n = 0;
	std::map<size_t,bool> used;

	if ( c  >  0 ) {
	  p_clusters_v = combined;
	  combined.clear();
	}

	
	for(unsigned i = 0; i < p_clusters_v.size(); ++i) {
	  used[i] = false;
	  
	  //is this necessary?
	  p_clusters_v[i].polygon()->UntanglePolygon();
	  
	}
      
	std::vector<size_t> overlaps; overlaps.reserve(p_clusters_v.size());

	for(unsigned c1 = 0; c1 < p_clusters_v.size(); ++c1) {

	  // std::cout << "{ ";

	  // for(unsigned i = 0; i < p_clusters_v.size(); ++i)
	  //   std::cout << "(" << i << "," << used[i] << "), ";

	  // std::cout << "}\n";
	
	  auto& p_cluster1 = p_clusters_v[c1];
	
	  overlaps.clear();
	
	  if ( used[c1] ) continue;
		  
	  for(unsigned c2 = 0; c2 < p_clusters_v.size(); ++c2) {
	  
	    if ( c1 == c2 )  continue;
	    if ( used[c2] )  continue;

	    auto& p_cluster2 = p_clusters_v[c2];
	  
	    //polygons overlap, get their points, convexHull on them, create new ProtoCluster
	    if ( p_cluster1.polygon()->PolyOverlapSegments( *( p_cluster2.polygon() ) ) ) { 
	      overlaps.push_back(c2);
	      n++;
	    }
	    
	  }


	  if(overlaps.size() == 0)
	    continue;
	
	  overlaps.push_back(c1);

	  int n_points = 0;
	  for(const auto& idx : overlaps)
	    n_points += p_clusters_v[idx].polygon()->Size();

	  //
	  std::vector<cv::Point> combine; combine.reserve ( n_points );
	  std::vector<cv::Point> hul;     hul.reserve     ( n_points );
	
	  std::vector<std::pair<float,float> > out;

	  for(const auto& idx : overlaps) {
	    for(unsigned p = 0; p < p_clusters_v[idx].polygon()->Size(); ++p) {
	      auto point = p_clusters_v[idx].polygon()->Point(p);
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
	
	  combined.emplace_back(out);

	  for(const auto& idx : overlaps)
	    used[idx] = true;

	} // end C1 loop
      
	for(const auto& u : used) {
	  if ( u.second ) continue;

	  std::vector<std::pair<float,float> > out ( p_clusters_v[u.first].polygon()->Size() );

	  //loop over points
	  for(unsigned k = 0; k < p_clusters_v[u.first].polygon()->Size(); ++k) 

	    out[k] = p_clusters_v[u.first].polygon()->Point(k);
	
	  combined.emplace_back(out);
	}

	if( n == 0 )
	  break;
	
	c++;
      }

      //that was terrible!
      
      std::cout << "\t==> Combined size...: " << combined.size() << "\n";

      //print out the hulls
      for(unsigned k = 0; k < combined.size(); ++k) {
	
	_hulls2.resize(combined.size());
	
	for(unsigned p = 0; p < combined[k].polygon()->Size(); ++p)  {

	  _hulls2[k].emplace_back(combined[k].polygon()->Point(p));

	}
      }
	
      //Ok now we can make "clusters"

      
      
	
          
    
  
      // Protoclusters are just convex hull objects in Polygon2D
      // They hold a pointer list of hough segments that begin, and pass through them
      // std::cout<< "\t(" << __FUNCTION__ << ")==> Assigning hough segments to protoclusters\n";
      // for(auto& proto : p_clusters_v) {
      // 	for(auto& hough : _houghs) {

      // 	  auto crossings = proto.polygon()->LineCross(std::make_pair(hough[0],hough[1]),
      // 						      std::make_pair(hough[2],hough[3]));
	  
      // 	  if( crossings > 0 )

      // 	    proto.AddLine(hough.data());
	  
      // 	}

      // 	//Determine most probable direction of proto cluster
      // 	//TODO
      // }
      
      
      //The ghetto way to see contours...
      _contour_v2.resize( cv_contour_v.size() );
      
      for(size_t c_index=0; c_index < cv_contour_v.size(); ++c_index) {
	
    	auto& cv_contour  = cv_contour_v[c_index];
	auto& out_contour = _contour_v2 [c_index];
	
    	out_contour.resize(cv_contour.size());
	
    	for(size_t p_index=0; p_index<cv_contour.size(); ++p_index)
	  out_contour[p_index] = std::make_pair(cv_contour[p_index].x,
						cv_contour[p_index].y);
      
      }
    
    }
    
    return evt_clusters;
  }


  size_t AlgoClusterHough::NumContours() const {
    return _contour_v.size();
  }
  
  PyObject* AlgoClusterHough::GetContour(const size_t contour_index) {

    if(contour_index > _contour_v.size()){
      std::cout << "\t==X Invalid contour ID requested...\n";
      throw std::exception();
    }

    int* dim_data = new int[2];
    dim_data[0] = 2;
    dim_data[1] = (int)(_contour_v[contour_index].size());

    return PyArray_FromDimsAndData( 2, dim_data, NPY_DOUBLE, (char*) &( _contour_v[contour_index].raw_data()[0] ) );
    
  }

    

}
#endif
