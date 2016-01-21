/**
 * \file Polygon.h
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief 2D polygon object
 *
 * @author kazuhiro & david caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/

#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <utility>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>


#include "LArUtil/GeometryHelper.h"

//a polygon is a vector of std::pairs with first = x coordinate
//and second = y coordinate of that vertex
//access vertices with Point function. Points are:
//0   = first ordered vertex
//n-1 = last ordered vertex (n=size of polygon)
//n   = first vertex again
//>n  = invalid...return error message

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//                BEGIN POLYGON CLASS               //
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
class Polygon{
  
 public:

  /// default constructor
  Polygon() { }
  /// constructor starting from list of edges for polygon
  Polygon(const std::vector< std::pair<float,float> > &points) { vertices = points; }
  /// Create Intersection Polygon from 2 polygons
  Polygon(const Polygon &poly1, const Polygon &poly2); 
  /// return number of edges in polygon
  unsigned int Size() const { return vertices.size(); } 
  /// return pair with information for a specific edge
  const std::pair<float,float>& Point(unsigned int p) const; 
  /// Project Polygon applying translation and rotation
  std::pair<float,float> Project(const std::pair<float,float> &p,float theta) const;
  /// Return Polygon Area
  float Area() const;
  /// return polygon perimeter
  float Perimeter() const;
  /// boolean: do these polygons overlap?
  bool PolyOverlap(const Polygon &poly2) const;
  /// boolean: is a point inside the polygon?
  bool PointInside(const std::pair<float,float> &point) const;
  /// boolean: line crosses polygon?
  int  LineCross(const std::pair<float,float> &point1,
		 const std::pair<float,float> &point2) const;
  /// check if poly2 is fully contained in poly1
  bool Contained(const Polygon &poly2) const; 
  /// untangle polygon
  void UntanglePolygon();
  /// clear polygon's points
  void Clear() { vertices.clear(); }
  
  ///Calculate the opening angle at the specified vertex:
  float InteriorAngle(unsigned int p) const;

  /// Return minimum distance between two polygons
  float Distance(const Polygon &poly2);
  
  
  friend bool operator==(const Polygon& lhs, const Polygon& rhs);
  friend bool operator!=(const Polygon& lhs, const Polygon& rhs);
  friend std::ostream &operator<<(std::ostream &out, Polygon poly);     //output

  
  /// utility function used in polygon overlap determination
  bool PolyOverlapSegments(const Polygon &poly2) const;
  
  
private:

  /// vector listing the polygon edges
  std::vector< std::pair<float,float> > vertices;

  /// utility function used by PolyOverlap to determine overlap
  bool Overlap(float slope, const Polygon &poly2, const std::pair<float,float> &origin) const;



};


/** @} */ // end of doxygen group

#endif
