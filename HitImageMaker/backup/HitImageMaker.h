/**
 * \file HitImageMaker.h
 *
 * \ingroup Play
 * 
 * \brief Class def header for a class HitImageMaker
 *
 * @author kazuhiro
 */

/** \addtogroup Play

    @{*/

#ifndef LARLITE_HITIMAGEMAKER_H
#define LARLITE_HITIMAGEMAKER_H

#include "Analysis/ana_base.h"
#include "opencv2/core/mat.hpp"
#include "Utils/DataTypes.h"
struct _object;
typedef _object PyObject;

#ifndef __CINT__
#include "Python.h"
#include "numpy/arrayobject.h"
#endif

namespace larlite {


  class HitImageMaker : public ana_base{
  
  public:

    /// Default constructor
    HitImageMaker();

    /// Default destructor
    virtual ~HitImageMaker(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();


    PyObject* GetCanny(const size_t plane);
    PyObject* GetContour(const size_t plane, const size_t contour_index);
    size_t NumContours(const size_t plane) const;

    ::cv::Mat* GetMat(const size_t plane);
    
  protected:
    void init();
    std::vector<::cv::Mat> _mat_v;
    std::vector<::cv::Mat> _canny_v;
    
    std::vector< std::vector<larcv::Point2DArray> >_contour_v;

    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
