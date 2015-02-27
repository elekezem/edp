/**************************************************************************
 *   mathtools.h                                                          *
 *                                                                        *
 *   EDP                                                                  *
 *                                                                        *
 *   Authors: Ivo Filot                                                   *
 *            Bart Zijlstra                                               *
 *            Emiel Hensen                                                *
 *                                                                        *
 *   (C) Copyright 2015 Inorganic Materials Chemistry                     *
 *                                                                        *
 *   This is a legal licensing agreement (Agreement) between              *
 *   You (an individual or single legal entity) and                       *
 *   Inorganic Materials Chemistry (IMC) governing the in-house use       *
 *   of the EDP software product (Software).                              *
 *   By downloading, installing, or using Software, You agree to be bound *
 *   by the license terms as given in the LICENSE file                    *
 *                                                                        *
 **************************************************************************/

/*
 * Vector class in 3D space
 */
class Vector {
private:
  float r[3];
public:
  Vector();
  Vector(float x, float y, float z);
private:
};

/*
 * 3x3 Matrix class
 */
class Matrix {
private:
  float* r[3];
public:
  Matrix();
  ~Matrix();
  void inverse();
  float* operator[](const unsigned int &i);
  const float* operator[](const unsigned int &i) const;
private:
};

/*
 * Plane class
 */
class Plane {
private:
  Vector origin; // origin of the plane
  Vector normal; // normal vector
public:
  Plane(const Vector &r, const Vector &n);
private:

};