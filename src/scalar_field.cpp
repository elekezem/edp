/**************************************************************************
 *   vasp_parser.cpp                                                      *
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

#include "scalar_field.h"

ScalarField::ScalarField(const std::string &_filename) {
  this->filename = _filename;
  this->scalar = -1;
}

void ScalarField::copy_grid_dimensions(unsigned int _grid_dimensions[]) {
  for(unsigned int i=0; i<3; i++) {
    _grid_dimensions[i] = this->grid_dimensions[i];
  }
}

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

XYZ ScalarField::realspace_to_grid(const double &i,
                     const double &j,
                     const double &k) const {
  XYZ r;
  r.x = imat[0][0] * i + imat[0][1] * j + imat[0][2] * k;
  r.y = imat[1][0] * i + imat[1][1] * j + imat[1][2] * k;
  r.z = imat[2][0] * i + imat[2][1] * j + imat[2][2] * k;

  r.x *= float(this->grid_dimensions[0]);
  r.y *= float(this->grid_dimensions[1]);
  r.z *= float(this->grid_dimensions[2]);

  return r;
}

/* read the scalar value at the 2nd line of the chgcar file */
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

/* read the unit cell matrix from the chgcar file */
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

/* read the list of atoms from the chgcar file */
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

void ScalarField::read_grid_dimensions(bool debug) {
  if(debug) std::cout << "Reading grid dimensions...\t\t";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  for(unsigned int i=0; i<7 + this->nrat.size() + 2; i++) {
    std::getline(infile, line);
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
}

void ScalarField::read_grid(bool debug) {
  std::cout.setf(std::ios_base::unitbuf); // flush after every "<<"
  if(debug) std::cout << "Reading grid values...";
  std::ifstream infile(this->filename.c_str());
  std::string line;
  for(unsigned int i=0; i<7 + this->nrat.size() + 2; i++) {
    std::getline(infile, line);
  }

  this->gridsize = this->grid_dimensions[0] * this->grid_dimensions[1]
                         * this->grid_dimensions[2];
  this->gridptr = new float[this->gridsize];  // spin up
  this->gridptr2 = new float[this->gridsize]; // spin down

  /* read spin up */
  unsigned int i=0;
  unsigned int cur=0; // for the counter
  unsigned int linecounter=0; // for the counter
  unsigned int wordcounter=0; // for the counter
  pcrecpp::RE re("([0-9E.+-]+)");
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
    while(re.FindAndConsume(&input, &this->gridptr[i])) {
      i++;
      wordcounter++;
    }
    if(wordcounter < 5) {
      std::cout << "I think I am not reading a VASP CHGCAR file, I am stopping." << std::endl;
      break;
    }


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

/* output everything to std::cout */
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

float ScalarField::get_value_interp(const float &x, const float &y, const float &z) {
  if(x > this->get_max_direction(0) || y > this->get_max_direction(1) || z > this->get_max_direction(2)) {
    return 0.0;
  }
  if(x < 0 || y < 0 || z < 0) {
    return 0.0;
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

  float c00 = this->get_value(x0, y0, z0) * (1-xd) +
              this->get_value(x1, y0, z0) * (1-xd);
  float c10 = this->get_value(x0, y1, z0) * (1-xd) +
              this->get_value(x1, y1, z0) * (1-xd);
  float c01 = this->get_value(x0, y0, z1) * (1-xd) +
              this->get_value(x1, y0, z1) * (1-xd);
  float c11 = this->get_value(x0, y1, z1) * (1-xd) +
              this->get_value(x1, y1, z1) * (1-xd);

  float c0 = c00 * (1 - yd) + c10 * yd;
  float c1 = c01 * (1 - yd) + c11 * yd;

  float c = c0 * (1 - zd) + c1 * zd;

  return c;
}

float ScalarField::get_max_direction(const unsigned int &dim) {
  float sum = 0;
  for(unsigned int i=0; i<3; i++) {
    sum += this->mat[i][dim];
  }
  return sum;
}

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

ScalarField::~ScalarField() {
  delete[] this->gridptr;
  delete[] this->gridptr2;
}