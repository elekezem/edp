/**************************************************************************
 *   edp.cpp                                                              *
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

#include <iostream>
#include "mathtools.h"

int main() {
  std::cout << "Running EDP" << std::endl;

  // assuming a square box of 10x10x10 A, construct a plane given by
  Vector r(5.0, 5.0, 5.0);
  // and a direction
  Vector n(0.0, 1.0, 0.0);
  // that constructs a plane
  Plane p(r, n);

  return 0;
}