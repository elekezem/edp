/**************************************************************************
 *   unitcell.cpp                                                         *
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

#include "unitcell.h"

Unitcell::Unitcell(const Vector &x, const Vector &y, const Vector &z) {
	for(unsigned int i=0; i<3; i++) {
		this->m[0][i] = x[i];
		this->m[1][i] = y[i];
		this->m[2][i] = z[i];
	}
}

float* Unitcell::operator[](const unsigned int &i) {
	return this->m[i];
}
 
const float* Unitcell::operator[](const unsigned int &i) const {
	return this->m[i];
}

void plot_unitcell_unfold(Plotter &plt, const Unitcell &u) {
	Color black(0,0,0);

	plt.draw_line(500, 500, 500 + u[0][0] * 10, 500 - u[0][1] * 10, black, 1.0);
	plt.draw_line(500, 500, 500 + u[1][0] * 10, 500 - u[1][1] * 10, black, 1.0);
	plt.draw_line(500 + u[0][0] * 10, 500 - u[0][1] * 10, 500 + (u[0][0] + u[1][0]) * 10, 500 - (u[0][1] + u[1][1]) * 10, black, 1.0);
	plt.draw_line(500 + u[1][0] * 10, 500 - u[1][1] * 10, 500 + (u[0][0] + u[1][0]) * 10, 500 - (u[0][1] + u[1][1]) * 10, black, 1.0);
}