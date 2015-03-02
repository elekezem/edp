/**************************************************************************
 *   unitcell.h                                                           *
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

#ifndef _UNITCELL_H
#define _UNITCELL_H

#include "mathtools.h"
#include "plotter.h"

class Unitcell {
private:
  Matrix m;
public:
  Unitcell(const Vector &x, const Vector &y, const Vector &z);
  float* operator[](const unsigned int &i);
  const float* operator[](const unsigned int &i) const;
private:
};

void plot_unitcell_unfold(Plotter &plt, const Unitcell &u);

#endif //_UNITCELL_H