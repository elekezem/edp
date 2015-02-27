#include "mathtools.h"

Vector::Vector() {
  this->r[0] = 0;
  this->r[1] = 0;
  this->r[2] = 0;
}

Vector::Vector(float x, float y, float z) {
  this->r[0] = x;
  this->r[1] = y;
  this->r[2] = z;
}

Matrix::Matrix() {
  for(unsigned int i=0; i<3; i++) {
    this->r[i] = new float[3];
    for(unsigned int j=0; j<3; j++) {
      this->r[i][j] = 0.0f;
    }
  }

}

Matrix::~Matrix() {
  for(unsigned int i=0; i<3; i++) {
    delete[] this->r[i];
  }
}