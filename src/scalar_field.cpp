/**************************************************************************
 *   scalar_field.cpp                                                     *
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

#include "scalar_field.h"

/*
 * Default constructor
 *
 * Usage: ScalarField sf("CHGCAR");
 *
 * Constructs a ScalarField by pointing it at a file on the HD.
 * The existence of the file is *not* being checked
 */
ScalarField::ScalarField(const std::string &_filename) {
  this->filename = _filename;
  this->scalar = -1;
}

/*
 * Destructor
 *
 * Empties the arrays
 *
 */
ScalarField::~ScalarField() {
  delete[] this->gridptr;
  delete[] this->gridptr2;
}

/*
 * void output()
 *
 * Outputs a summary of the ScalarField to std::cout.
 * Mainly used for debugging purposes.
 *
 */
void ScalarField::output() const {
  std::cout << "Scalar: "<< this->scalar << std::endl;
  std::cout << std::endl;
  std::cout << "Matrix: ";
  for(unsigned i=0; i<3; i++) {
    for(unsigned j=0; j<3; j++) {
      std::cout << this->mat[i][j] << "\t";
    }
    std::cout << std::endl;
    std:: cout << "\t";
  }
  std::cout << std::endl;
  std::cout << "Inverse: ";
  for(unsigned i=0; i<3; i++) {
    for(unsigned j=0; j<3; j++) {
      std::cout << this->imat[i][j] << "\t";
    }
    std::cout << std::endl;
    std:: cout << "\t";
  }
  std::cout << std::endl;
  std::cout << "ion types: " << this->nrat.size() << " ( ";
  for(unsigned i=0; i<this->nrat.size(); i++) {
    std::cout << this->nrat[i] << " ";
  }
  std::cout << ")" << std::endl;
  std::cout << std::endl;
  std::cout << "Grid dimensions: ";
  for(unsigned i=0; i<3; i++) {
    std::cout << this->grid_dimensions[i] << "\t";
  }
  std::cout << std::endl;
}

/*
 * void read(bool debug)
 *
 * Wrapper function that reads in the OUTCAR file
 *
 * Usage: sf.read(true);
 *
 * Currently, only VASP4 output is being supported!
 */
void ScalarField::read(bool debug) {
  this->read_scalar(debug);
  this->read_matrix(debug);
  this->read_atoms(debug);
  this->read_grid_dimensions(debug);
  this->read_grid(debug);
}

/*
 * void read_scalar(bool debug)
 *
 * Read the scalar value from the 2nd line of the
 * CHGCAR file. Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_scalar(bool debug) {
  if(debug) std::cout << "Reading scalar...\t\t\t";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  std::getline(infile, line); // discard this line

  std::getline(infile, line);
  double val = -1;
  pcrecpp::RE re("^\\s*([0-9.-]+)\\s*$");
  re.FullMatch(line , &val);
  this->scalar = val;
  if(debug) std::cout << "[Done]" << std::endl;
}

/*
 * void read_matrix(bool debug)
 *
 * Reads the matrix that defines the unit cell
 * in the CHGCAR file. The inverse of that matrix
 * is automatically constructed.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_matrix(bool debug) {
  if(debug) std::cout << "Reading unitcell matrix...\t\t";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  for(unsigned int i=0; i<2; i++) { // discard first two lines
    std::getline(infile, line);
  }

  // setup match pattern
  pcrecpp::RE re("^\\s*([0-9.-]+)\\s+([0-9.-]+)\\s+([0-9.-]+)\\s*$");
  for(unsigned int i=0; i<3; i++) {
    std::getline(infile, line);
    re.FullMatch(line , &mat[i][0], &mat[i][1], &mat[i][2]);
  }

  for(unsigned int i=0; i<3; i++) {
    for(unsigned int j=0; j<3; j++) {
      this->mat[i][j] *= this->scalar;
    }
  }

  // also construct inverse matrix
  this->calculate_inverse();

  if(debug) std::cout << "[Done]" << std::endl;
}

/*
 * void read_atoms(bool debug)
 *
 * Read the number of atoms of each element. These
 * numbers are used to skip the required amount of
 * lines.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_atoms(bool debug) {
  if(debug) std::cout << "Reading atoms...\t\t\t";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  for(unsigned int i=0; i<5; i++) { // discard first two lines
    std::getline(infile, line);
  }

  int val = 0;
  std::getline(infile, line);
  pcrecpp::RE re("([0-9]+)");
  pcrecpp::StringPiece input(line);
  while(re.FindAndConsume(&input, &val)) {
    this->nrat.push_back(val);
  }
  if(debug) std::cout << "[Done]" << std::endl;
}

/*
 * void read_grid_dimensions(bool debug)
 *
 * Read the number of gridpoints in each
 * direction.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_grid_dimensions(bool debug) {
  if(debug) std::cout << "Reading grid dimensions...\t\t";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  // skip lines that contain atoms
  for(unsigned int i=0; i<9; i++) {
    std::getline(infile, line);
  }
  for(unsigned int i=0; i<this->nrat.size(); i++) {
    for(unsigned int j=0; j<this->nrat[i]; j++) {
        std::getline(infile, line);
    }
  }

  this->gridline = line;

  pcrecpp::RE re("([0-9]+)");
  pcrecpp::StringPiece input(line);
  unsigned int i = 0;
  unsigned int val = 0;
  while(re.FindAndConsume(&input, &val)) {
    this->grid_dimensions[i] = val;
    i++;
  }
  if(debug) std::cout << "[Done]" << std::endl;
  if(debug) std::cout << "GRID: " << this->grid_dimensions[0] << "x" <<
                                     this->grid_dimensions[1] << "x" <<
                                     this->grid_dimensions[2] << std::endl;
}

/*
 * void read_grid(bool debug)
 *
 * Read all the grid points. This function depends
 * on the the gridsize being set via the
 * read_grid_dimensions() function.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
void ScalarField::read_grid(bool debug) {
  std::cout.setf(std::ios_base::unitbuf); // flush after every "<<"
  if(debug) std::cout << "Reading grid values...";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  // skip lines that contain atoms
  for(unsigned int i=0; i<9; i++) {
    std::getline(infile, line);
  }
  for(unsigned int i=0; i<this->nrat.size(); i++) {
    for(unsigned int j=0; j<this->nrat[i]; j++) {
        std::getline(infile, line);
    }
  }

  this->gridsize = this->grid_dimensions[0] * this->grid_dimensions[1]
                         * this->grid_dimensions[2];
  this->gridptr = new float[this->gridsize];  // spin up
  this->gridptr2 = new float[this->gridsize]; // spin down (not being used now)

  /* read spin up */
  unsigned int i=0;
  unsigned int cur=0;         // for the counter
  unsigned int linecounter=0; // for the counter
  unsigned int wordcounter=0; // for the counter
  pcrecpp::RE re("([0-9Ee.+-]+)");
  pcrecpp::RE aug("augmentation.*");
  while(std::getline(infile, line)) {
    // stop looping when a second gridline appears (this
    // is where the spin down part starts)
    if(line.compare(this->gridline) == 0) {
      std::cout << "I am breaking the loop" << std::endl;
      break;
    }

    if(aug.FullMatch(line)) {
      break;
    }

    pcrecpp::StringPiece input(line);

    wordcounter = 0;
    if(i > this->gridsize) {
        break;
    }
    while(re.FindAndConsume(&input, &this->gridptr[i])) {
      i++;
      wordcounter++;
    }

    /*
     * Track the progress of the read procedure. (this is the task that takes the
     * most amount of time)
     */
    if(debug) {
      if(i > (this->gridsize / 5 * cur)) {
        std::cout << "[" << int(cur * 20) << " %]";
        cur++;
      }
    }
    linecounter++;
  }
  if(debug) std::cout << "[100%]" << std::endl;
  if(debug) std::cout << "End reading " << linecounter << " lines" << std::endl;
  if(debug) std::cout << "Grabbed " << i << "/" << this->gridsize << " values" << std::endl;

  // /* read spin down */
  // i=0;
  // while(std::getline(infile, line)) {
  //   pcrecpp::StringPiece input(line);
  //   while(re.FindAndConsume(&input, &this->gridptr2[i])) {
  //     i++;
  //   }
  // }
  // if(debug) std::cout << "[Done]" << std::endl;
}

/*
 * float get_value_interp(x,y,z)
 *
 * Grabs a value from the 3D scalar field. Calculate the value
 * by using a trilinear interpolation.
 *
 * The trilinear interpolation algorithm has been extracted from:
 * http://paulbourke.net/miscellaneous/interpolation/
 *
 * Future algorithm can make use of a cubic interpolation.
 *
 */
float ScalarField::get_value_interp(const float &x, const float &y, const float &z) {
  if(x > this->get_max_direction(0) || y > this->get_max_direction(1) || z > this->get_max_direction(2)) {
    return -1.0;
  }
  if(x < 0 || y < 0 || z < 0) {
    return -1.0;
  }

  // cast the input to the nearest grid point
  XYZ r = this->realspace_to_grid(x,y,z);

  // calculate value using trilinear interpolation
  float xd = remainderf(r.x, 1.0);
  float yd = remainderf(r.y, 1.0);
  float zd = remainderf(r.z, 1.0);

  if(xd < 0) xd += 1.0;
  if(yd < 0) yd += 1.0;
  if(zd < 0) zd += 1.0;

  float x0 = floor(r.x);
  float x1 = ceil(r.x);
  float y0 = floor(r.y);
  float y1 = ceil(r.y);
  float z0 = floor(r.z);
  float z1 = ceil(r.z);

  return
  this->get_value(x0, y0, z0) * (1.0 - xd) * (1.0 - yd) * (1.0 - zd) +
  this->get_value(x1, y0, z0) * xd         * (1.0 - yd) * (1.0 - zd) +
  this->get_value(x0, y1, z0) * (1.0 - xd) * yd         * (1.0 - zd) +
  this->get_value(x0, y0, z1) * (1.0 - xd) * (1.0 - yd) * zd         +
  this->get_value(x1, y0, z1) * xd         * (1.0 - yd) * zd         +
  this->get_value(x0, y1, z1) * (1.0 - xd) * yd         * zd         +
  this->get_value(x1, y1, z0) * xd         * yd         * (1.0 - zd) +
  this->get_value(x1, y1, z1) * xd         * yd         * zd;
}

/*
 * float get_max_direction(dim)
 *
 * Get the maximum value in a particular dimension. This is a convenience
 * function for the get_value_interp() function.
 *
 */
float ScalarField::get_max_direction(const unsigned int &dim) {
  float sum = 0;
  for(unsigned int i=0; i<3; i++) {
    sum += this->mat[i][dim];
  }
  return sum;
}

/*
 * void calculate_inverse()
 *
 * Calculates the inverse of a 3x3 matrix. This is a convenience
 * function for the read_matrix() function.
 *
 */
void ScalarField::calculate_inverse() {
  float det = 0;
  for(unsigned int i=0;i<3;i++) {
    det += (this->mat[0][i]*(this->mat[1][(i+1)%3]*this->mat[2][(i+2)%3] - this->mat[1][(i+2)%3]*this->mat[2][(i+1)%3]));
  }

  for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++) {
           this->imat[i][j] = ((this->mat[(i+1)%3][(j+1)%3] * this->mat[(i+2)%3][(j+2)%3]) - (this->mat[(i+1)%3][(j+2)%3]*this->mat[(i+2)%3][(j+1)%3]))/ det;
      }
   }
}

/*
 * float get_value(i,j,k)
 *
 * Grabs the value at a particular grid point.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
const float& ScalarField::get_value(const unsigned int i,
                       const unsigned int j,
                       const unsigned int k) const {
  if(i >= this->grid_dimensions[0]) {
    std::cout << "ERROR: Cannot access x=" << i << std::endl;
  }
  if(j >= this->grid_dimensions[1]) {
    std::cout << "ERROR: Cannot access y=" << j << std::endl;
  }
  if(k >= this->grid_dimensions[2]) {
    std::cout << "ERROR: Cannot access z=" << k << std::endl;
  }
  unsigned int idx = k * this->grid_dimensions[0] * this->grid_dimensions[1] +
                     j * this->grid_dimensions[0] +
                     i;
  if(idx > this->gridsize) {
    std::cout << "Trying to allocate value outside gridspace: (" << i << "," << j << "," << k << ")" << std::endl;
  }
  return this->gridptr[idx];
}

/*
 * XYZ grid_to_realspace(i,j,k)
 *
 * Converts a grid point to a realspace vector. This function
 * is not being used at the moment.
 *
 */
XYZ ScalarField::grid_to_realspace(const double &i,
                     const double &j,
                     const double &k) const {
  float dx = (double)i / (double)grid_dimensions[0];
  float dy = (double)j / (double)grid_dimensions[1];
  float dz = (double)k / (double)grid_dimensions[2];

  XYZ r;
  r.x = mat[0][0] * dx + mat[0][1] * dy + mat[0][2] * dz;
  r.y = mat[1][0] * dx + mat[1][1] * dy + mat[1][2] * dz;
  r.z = mat[2][0] * dx + mat[2][1] * dy + mat[2][2] * dz;

  return r;
}

/*
 * XYZ realspace_to_grid(i,j,k)
 *
 * Convert 3d realspace vector to a position on the grid. Non-integer
 * values (i.e. floating point) are given as the result.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
XYZ ScalarField::realspace_to_grid(const double &i,
                     const double &j,
                     const double &k) const {
  XYZ r;
  r.x = imat[0][0] * i + imat[0][1] * j + imat[0][2] * k;
  r.y = imat[1][0] * i + imat[1][1] * j + imat[1][2] * k;
  r.z = imat[2][0] * i + imat[2][1] * j + imat[2][2] * k;

  r.x *= float(this->grid_dimensions[0]-1);
  r.y *= float(this->grid_dimensions[1]-1);
  r.z *= float(this->grid_dimensions[2]-1);

  return r;
}
