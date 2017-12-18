//
//  implicit.cpp
//  cs175-final
//
//  Created by Gabe Montague on 12/16/17.
//  Copyright © 2017 cs175. All rights reserved.
//


#include "cvec.h"
#include "implicit.hpp"

CubeImplicit::CubeImplicit() {
}
  
CubeImplicit::CubeImplicit(Cvec3 center, double size, double thick) {
  m_center = center;
  m_size = size;
}
  
double CubeImplicit::getValue(const Cvec3 pos) const {
  const Cvec3 posRel = pos - m_center;
  const double absX = std::abs(posRel[0]);
  const double absY = std::abs(posRel[1]);
  const double absZ = std::abs(posRel[2]);
  
  if (absX > absY && absX > absZ) {
    return absX - m_size;
  } else if (absY > absZ) {
    return absY - m_size;
  } else {
    return absZ - m_size;
  }
}
  
Cvec3 CubeImplicit::getNorm(const Cvec3 pos) const {
  const Cvec3 posRel = pos - m_center;
  const double absX = std::abs(posRel[0]);
  const double absY = std::abs(posRel[1]);
  const double absZ = std::abs(posRel[2]);
  
  if (absX > absY && absX > absZ) {
    return Cvec3(posRel[0] > 0 ? 1 : -1, 0, 0);
  } else if (absY > absZ) {
    return Cvec3(0, posRel[1] > 0 ? 1 : -1, 0);
  } else {
    return Cvec3(0, 0, posRel[2] > 0 ? 1 : -1);
  }
}

SphereImplicit::SphereImplicit() {
}

SphereImplicit::SphereImplicit(Cvec3 center, double radius, double thick) {
  m_center = center;
  m_radius = radius;
}

double SphereImplicit::getValue(const Cvec3 pos) const {
  const Cvec3 posRel = pos - m_center;
  return norm(posRel) - m_radius;
}

Cvec3 SphereImplicit::getNorm(const Cvec3 pos) const {
  const Cvec3 posRelNorm = (pos - m_center).normalize();
  //return Cvec3f(posRelNorm[0], posRelNorm[1], posRelNorm[2]);
  return posRelNorm;
}

double TorusImplicit::getValue(const Cvec3 pos) const {
  const Cvec3 posRel = pos - m_center;
  Cvec3 posRelFlat = posRel;
  posRelFlat[1] = 0.0;
  double toSqr = m_radius - norm(posRelFlat);
  return toSqr * toSqr + posRel[1] * posRel[1] - m_thickness;
}

Cvec3 TorusImplicit::getNorm(const Cvec3 pos) const {
  const Cvec3 posRel = pos - m_center;
  Cvec3 posRelFlat = posRel;
  posRelFlat[1] = 0.0;
  
  const double normFlat = norm(posRelFlat);
  
  const double dx = posRel[0] * (2 - 2 * m_radius / normFlat);
  const double dy = 2 * posRel[1];
  const double dz = posRel[2] * (2 - 2 * m_radius / normFlat);
  
  const Cvec3 result = Cvec3(dx, dy, dz);
  const double normResult = norm(result);
  if (normResult < CS175_EPS) {
    return Cvec3();
  } else {
    return result / normResult;
  }
}

  
double ToggleImplicit::getValue(const Cvec3 pos) const {
  
  if (m_mode == 0) {
    return m_sphere.getValue(pos);
  } else if (m_mode == 2) {
    return m_cube.getValue(pos);
  } else if (m_mode == 1) {
    return m_torus.getValue(pos);
  } else {
    return m_intersection.getValue(pos);
  }
}
Cvec3 ToggleImplicit::getNorm(const Cvec3 pos) const {
  
  if (m_mode == 0) {
    return m_sphere.getNorm(pos);
  } else if (m_mode == 2) {
    return m_cube.getNorm(pos);
  } else if (m_mode == 1) {
    return m_torus.getNorm(pos);
  } else {
    return m_intersection.getNorm(pos);
  }
}


