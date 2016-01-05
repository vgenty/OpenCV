/**
 * \file ContourTreeMaker.h
 *
 * \ingroup Play
 * 
 * \brief Class def header for a class ContourTreeMaker
 *
 * @author kazuhiro
 */

/** \addtogroup Play

    @{*/
#ifndef CONTOURTREEMAKER_H
#define CONTOURTREEMAKER_H

#include <iostream>
#include <TTree.h> // For TTree
#include <vector>  // For std::vector
#include <utility> // For std::pair

namespace imutil {

  typedef std::pair<double,double>          Coordinate_t;
  typedef std::vector<imutil::Coordinate_t> Contour_t;
  typedef std::vector<imutil::Contour_t>    ContourArray_t;
  typedef size_t ContourID_t;
  
  /**
     \class ContourTreeMaker
     Simple class to create & fill TTree
  */
  class ContourTreeMaker{
    
  public:
    
    /// Default constructor
    ContourTreeMaker(const std::string tree_name = "tree" );
    
    /// Default destructor
    ~ContourTreeMaker(){}

    /// Contour-fill method 1: just get contour reference and fill by yourself.
    imutil::ContourArray_t& ContourArray() { return _contour_array; }

    /// Contour-fill method 2: fill via function
    void FillContour(const ContourID_t contour_index, const double x, const double y);

    /// TTree-fill method
    void Fill() { _tree.Fill(); }

    /// TTree getter
    const TTree& Tree() { return _tree; }
    
  protected:

    imutil::ContourArray_t _contour_array;
    TTree _tree;

  };
}

#endif
/** @} */ // end of doxygen group 

