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
#include "scalar_field.h"

int main() {
  std::cout << "Running EDP" << std::endl;

  ScalarField sf("CHGCAR");
  sf.read_scalar(true);
  sf.read_matrix(true);
  sf.read_atoms(true);
  sf.read_grid_dimensions(true);
  sf.read_grid(true);
  sf.output();
//sf.get_value_interp(5.021,5.021,5.021);

  // create a canvas
  ColorScheme scheme(-10,10);
  Plotter plt(1000, 1000);
  for(float x=0; x<10; x+=0.01) {
    for(float y=0; y<10; y+=0.01) {
      plt.draw_filled_rectangle(x * 100, y * 100, 1, 1, scheme.get_color(sf.get_value_interp(x,5,y), true));
    }
  }
  plt.write("test.png");

  // // assuming a unit cell of 10x10x10
  // Unitcell u(Vector(10,0,0), Vector(0,10,0), Vector(0,0,10));

  // Vector r(5.0, 5.0, 5.0); 	// construction position
  // Vector n(0.0, 1.0, 0.0);	// construction normal
  // Plane p(r, n);

  // plot_unitcell_unfold(plt, u);

  // plt.write("test.png");

  return 0;
}