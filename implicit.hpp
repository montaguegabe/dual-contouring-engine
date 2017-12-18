//
//  implicit.h
//  cs175-final
//
//  Created by Gabe Montague on 12/16/17.
//  Copyright © 2017 cs175. All rights reserved.
//

#ifndef implicit_h
#define implicit_h

#include "cvec.h"

// All have same methods

class CubeImplicit {
public:
  
  CubeImplicit();
  
  CubeImplicit(Cvec3 center, double size, double thick);
  
  Cvec3 m_center;
  double m_size;
  
  double getValue(const Cvec3 pos) const;
  Cvec3 getNorm(const Cvec3 pos) const;
  
};

class SphereImplicit {
public:
  
  SphereImplicit();
  
  SphereImplicit(Cvec3 center, double radius, double thick);
  
  Cvec3 m_center;
  double m_radius;
  
  double getValue(const Cvec3 pos) const;
  Cvec3 getNorm(const Cvec3 pos) const;
  
};

class TorusImplicit {
public:
  
  TorusImplicit() {};
  
  TorusImplicit(Cvec3 center, double radius, double thickness) {
    m_center = center;
    m_radius = radius;
    m_thickness = thickness;
  }
  
  Cvec3 m_center;
  double m_radius;
  double m_thickness;
  
  double getValue(const Cvec3 pos) const;
  Cvec3 getNorm(const Cvec3 pos) const;
  
};

template <typename A, typename B>
class IntersectionImplicit {
public:
  IntersectionImplicit(Cvec3 center, double radius, double thick) : m_aField(center, radius, thick), m_bField(center, radius, thick) {
    
  }
  A m_aField;
  B m_bField;
  double getValue(const Cvec3 pos) const {
    double aVal = m_aField.getValue(pos);
    double bVal = -m_bField.getValue(pos);
    return aVal > bVal ? aVal : bVal;
  }
  Cvec3 getNorm(const Cvec3 pos) const {
    double aVal = m_aField.getValue(pos);
    double bVal = -m_bField.getValue(pos);
    /*if (std::abs(aVal - bVal) < 0.01) {
      return normalize(m_aField.getNorm(pos) - m_bField.getNorm(pos));
    }*/
    return aVal > bVal ? m_aField.getNorm(pos) : -m_bField.getNorm(pos);
  }
};

class ToggleImplicit {
public:
  ToggleImplicit(Cvec3 center, double radius, double thick, int shape) :
    m_sphere(center, radius, thick),
    m_cube(center, radius, thick),
    m_torus(center, radius, thick),
    m_intersection(center, radius, thick),
    m_mode(shape){
  }
  
  SphereImplicit m_sphere;
  CubeImplicit m_cube;
  TorusImplicit m_torus;
  IntersectionImplicit<SphereImplicit, TorusImplicit> m_intersection;
  int m_mode = 0;
  
  double getValue(const Cvec3 pos) const;
  Cvec3 getNorm(const Cvec3 pos) const;
  
};

#endif /* implicit_h */
