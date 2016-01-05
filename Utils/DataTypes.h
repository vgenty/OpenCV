#ifndef __OPENCV_UTILS_DATATYPES_H__
#define __OPENCV_UTILS_DATATYPES_H__

#include <vector>

namespace larcv {

  struct Point2D {
    double x, y;
    Point2D(double xv=0, double yv=0) : x(xv), y(yv) {}
  };

  class Point2DArray {

  public:
    Point2DArray() : _data() {}
    ~Point2DArray() {}
    
    void resize    (const size_t i) { _data.resize(i*2);  }
    void reserve   (const size_t i) { _data.reserve(i*2); }
    void push_back (const double x, const double y)
    { _data.push_back(x); _data.push_back(y); }
    void set       (const size_t i, const double x, const double y)
    { _data[i] = x; _data[_data.size()/2+i] = y; }

    size_t size   ()               const { return _data.size() / 2;               }
    double x      (const size_t i) const { return _data[i];                       }
    double y      (const size_t i) const { return _data[_data.size()/2+i];        }
    Point2D point (const size_t i) const
    { return Point2D(_data[i], _data[_data.size()/2+i]); }

    std::vector<double>& raw_data() { return _data; }
    //const std::vector<double>& raw_data() const { return _data; }
    
  private:
    std::vector<double> _data;

  };

}

#endif
