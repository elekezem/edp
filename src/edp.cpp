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
#include "unitcell.h"

int main() {
  std::cout << "Running EDP" << std::endl;

  // create a canvas
  Plotter plt(1000, 1000);

  // assuming a unit cell of 10x10x10
  Unitcell u(Vector(10,0,0), Vector(0,10,0), Vector(0,0,10));

  Vector r(5.0, 5.0, 5.0); 	// construction position
  Vector n(0.0, 1.0, 0.0);	// construction normal
  Plane p(r, n);

  plot_unitcell_unfold(plt, u);

  plt.write("test.png");

  return 0;
}