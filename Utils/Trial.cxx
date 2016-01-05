#ifndef TRIAL_CXX
#define TRIAL_CXX

#include "Trial.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "NDArrayConverter.h"

PyObject* Trial::GetImage() const
{
  cv::Mat image;
  image = cv::imread("/Users/kazuhiro/Downloads/Lenna.png", CV_LOAD_IMAGE_COLOR);   // Read the file
  
  if(! image.data )                              // Check for invalid input
    {
      std::cout <<  "Could not open or find the image" << std::endl ;
      return nullptr;
    }

  ::larcv::convert::NDArrayConverter converter;
  return converter.toNDArray(image);
}
#endif
