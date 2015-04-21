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

/*
 * The Electron Density Plotter (EDP) projects the 3D electron density onto
 * a 2D plane (i.e. a surface cut).
 *
 * The a program reads one or several CHGCAR files, performs elementwise
 * some mathematical operations on the content and creates a memory object
 * of the 3D scalar field.
 *
 * From this 3D scalar field, a surface cut is produced using a trilinear
 * interpolation routine.
 *
 */

#include <iostream>
#include "mathtools.h"
#include "scalar_field.h"
#include "planeprojector.h"

int main() {
    std::cout << "Running EDP" << std::endl;

    // read in field
    ScalarField sf("CHGCAR");
    sf.read(true);

    // define vectors and start points
    // vectors have to be normalized!
    Vector v1(1,0,0);
    Vector v2(0,0,1);
    Vector s(5,5,5);
    float scale = 200;

    // define intervals in Angstrom
    float li = -5.0;
    float hi = 5.0;

    float lj = -5.0;
    float hj = 5.0;

    PlaneProjector pp(&sf, 0, 5);
    pp.extract(v1, v2, s, scale, li, hi, lj, hj);
    pp.plot();
    pp.isolines(5);
    pp.write("test.png");

    return 0;
}
