#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>
#include <array>
#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS175_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r) {
    t_ = t;
    r_ = r;
  }

  explicit RigTForm(const Cvec3& t) {
    t_ = t;
    r_ = Quat();
  }

  explicit RigTForm(const Quat& r) {
    t_ = {0, 0, 0};
    r_ = r;
  }

  explicit RigTForm(const std::array<double, 7> posQuat) {
    t_ = {posQuat[0], posQuat[1], posQuat[2]};
    r_ = Quat(posQuat[3], posQuat[4], posQuat[5], posQuat[6]);
  }
  
  std::array<double, 7> getPosQuat() const {
    return {t_[0], t_[1], t_[2], r_.q_[0], r_.q_[1], r_.q_[2], r_.q_[3]};
  }
  
  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  Cvec4 operator * (const Cvec4& a) const {

    const Cvec3 inCoords = {a[0], a[1], a[2]};
    const bool translate = a[3];
    
    // Do rotation
    Cvec3 outCoords = r_ * inCoords;
    
    // Do translation
    if (translate) {
      outCoords += t_;
    }
    
    return {outCoords[0], outCoords[1], outCoords[2], translate};
  }

  RigTForm operator * (const RigTForm& a) const {
    const Quat outR = r_ * a.r_;
    const Cvec3 outT = t_ + r_ * a.t_;
    return RigTForm(outT, outR);
  }
};

inline RigTForm inv(const RigTForm& tform) {
  const Quat outR = inv(tform.getRotation());
  const Cvec3 outT = -(outR * tform.getTranslation());
  return RigTForm(outT, outR);
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
  Matrix4 t = Matrix4::makeTranslation(tform.getTranslation());
  Matrix4 r = quatToMatrix(tform.getRotation());
  return  t * r;
}

#endif
