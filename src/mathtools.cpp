/**************************************************************************
 *   mathtools.cpp                                                        *
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

#include "mathtools.h"

/*
 * Default vector constructor
 */
Vector::Vector() {
  this->r[0] = 0;
  this->r[1] = 0;
  this->r[2] = 0;
}

/*
 * Specific vector constructor
 */
Vector::Vector(float x, float y, float z) {
  this->r[0] = x;
  this->r[1] = y;
  this->r[2] = z;
}

/*
 * Allocator method to place a value in the Vector
 */
float Vector::operator[](const unsigned int &i) {
  return this->r[i];
}

/*
 * Allocator method to grab a value from the Vector
 */
const float& Vector::operator[](const unsigned int &i) const {
  return this->r[i];
}

float Vector::length() {
    return sqrt(this->r[0] * this->r[0] +
                this->r[1] * this->r[1] +
                this->r[2] * this->r[2]);
}

void Vector::normalize() {
    float length = this->length();
    this->r[0] /= length;
    this->r[1] /= length;
    this->r[2] /= length;
}

/*
 * Default matrix constructor
 */
Matrix::Matrix() {
  for(unsigned int i=0; i<3; i++) {
    this->r[i] = new float[3];
    for(unsigned int j=0; j<3; j++) {
      this->r[i][j] = 0.0f;
    }
  }
}

/*
 * Matrix destructor
 */
Matrix::~Matrix() {
  for(unsigned int i=0; i<3; i++) {
    delete[] this->r[i];
  }
}

/*
 * Allocator method to place a value in the matrix
 */
float* Matrix::operator[](const unsigned int &i) {
  return this->r[i];
}

/*
 * Allocator method to grab a value from the matrix
 */
const float* Matrix::operator[](const unsigned int &i) const {
  return this->r[i];
}

Plane::Plane(const Vector &r, const Vector &n) {
  this->origin = r;
  this->normal = n;
}

void Plane::parametrize() {

}
