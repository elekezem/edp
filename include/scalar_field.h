/**************************************************************************
 *   vasp_parser.h                                                        *
 *                                                                        *
 *   ISOTRON                                                              *
 *                                                                        *
 *   This program is free software; you can redistribute it and/or modify *
 *   it under the terms of the GNU General Public License as published by *
 *   the Free Software Foundation, version 2                              *
 *                                                                        *
 *   This program is distributed in the hope that it will be useful, but  *
 *   WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    *
 *   General Public License for more details.                             *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program; if not, write to the Free Software          *
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA            *
 *   02110-1301, USA.                                                     *
 *                                                                        *
 **************************************************************************/

#ifndef _SCALAR_FIELD_H
#define _SCALAR_FIELD_H

#include <string>
#include <vector>
#include <iostream>
#include <ios>
#include <sstream>
#include <fstream>
#include <pcrecpp.h>
#include <math.h> 
#include "xyz.h"

class ScalarField{
private:
	std::string filename;
	double scalar;
	double mat[3][3];
	double imat[3][3];

	unsigned int grid_dimensions[3];
	std::vector<unsigned int> nrat;
	std::string gridline;
	float* gridptr; // grid to first pos of float array
	float* gridptr2; // grid to first pos of float array
	unsigned int gridsize;

public:
	ScalarField(const std::string &_filename);
	void output() const;
	void copy_grid_dimensions(unsigned int _grid_dimensions[]);
	const float& get_value(const unsigned int i,
                       const unsigned int j,
                       const unsigned int k) const;
	XYZ grid_to_realspace(const double &i,
                       const double &j,
                       const double &k) const;
	XYZ realspace_to_grid(const double &i,
                       const double &j,
                       const double &k) const;
	~ScalarField();

public:
	void read_scalar(bool debug);
	void read_matrix(bool debug);
	void read_grid_dimensions(bool debug);
	void read_grid(bool debug);
	void read_atoms(bool debug);
	//void read_coordinates();

	float get_value_interp(const float &x, const float &y, const float &z);

private:
	float get_max_direction(const unsigned int &dim);
	void calculate_inverse();
};

#endif