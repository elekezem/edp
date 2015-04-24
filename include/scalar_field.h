/**************************************************************************
 *   scalar_field.h                                                       *
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
    double mat[3][3];   // matrix dimensions
    double imat[3][3];  // inverse of matrix

    unsigned int grid_dimensions[3];
    std::vector<unsigned int> nrat;
    std::string gridline;
    float* gridptr;  // grid to first pos of float array
    float* gridptr2; // grid to first pos of float array
    unsigned int gridsize;
    bool vasp5_input;

public:
    ScalarField(const std::string &_filename);
    void output() const;
    ~ScalarField();

public:
    void read(bool debug);

    /*
     * function for reading in the CHGCAR file
     */
private:
    void test_vasp5(bool debug);
    void read_scalar(bool debug);
    void read_matrix(bool debug);
    void read_grid_dimensions(bool debug);
    void read_grid(bool debug);
    void read_atoms(bool debug);

    /*
     * output and handler functions
     */
public:
    float get_value_interp(const float &x, const float &y, const float &z);

    /*
     * utility functions
     */
private:
    float get_max_direction(const unsigned int &dim);
    void calculate_inverse();

    /*
     * value extraction and dimensionality manipulators
     */
    const float& get_value(const unsigned int i,
                       const unsigned int j,
                       const unsigned int k) const;
    XYZ grid_to_realspace(const double &i,
                       const double &j,
                       const double &k) const;
    XYZ realspace_to_grid(const double &i,
                       const double &j,
                       const double &k) const;
    XYZ realspace_to_direct(const double &i,
                       const double &j,
                       const double &k) const;
};

#endif
