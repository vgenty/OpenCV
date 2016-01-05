/**
 * \file Trial.h
 *
 * \ingroup Utils
 * 
 * \brief Class def header for a class Trial
 *
 * @author kazuhiro
 */

/** \addtogroup Utils

    @{*/
#ifndef TRIAL_H
#define TRIAL_H

struct _object;
typedef _object PyObject;

#ifndef __CINT__
#include "Python.h"
#endif

/**
   \class Trial
   User defined class Trial ... these comments are used to generate
   doxygen documentation!
 */
class Trial{

public:

  /// Default constructor
  Trial(){}

  /// Default destructor
  ~Trial(){}

  PyObject* GetImage() const;

};

#endif
/** @} */ // end of doxygen group 

