/**************************************************************************
 *   planeprojector.h                                                     *
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

#ifndef _PLANEPROJECTOR_H
#define _PLANEPROJECTOR_H

#include <algorithm>
#include "plotter.h"
#include "mathtools.h"
#include "scalar_field.h"

class PlaneProjector {
private:
    ColorScheme* scheme;
    ScalarField* sf;
public:
    PlaneProjector(ScalarField* _sf);
    void plot(Vector _v1, Vector _v2, Vector _s, float _scale, float li, float hi, float lj, float hj);
private:
};

#endif
