/**************************************************************************
 *   planeprojector.cpp                                                   *
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

 #include "planeprojector.h"


PlaneProjector::PlaneProjector(ScalarField* _sf, float _min, float _max) {
    this->min = _min;
    this->max = _max;
    this->scheme = new ColorScheme(_min,_max);
    this->sf = _sf;
}

void PlaneProjector::extract(Vector _v1, Vector _v2, Vector _s, float _scale, float li, float hi, float lj, float hj, bool negative_values) {

    this->ix = int((hi - li) * _scale);
    this->iy = int((hj - lj) * _scale);

    std::cout << "Creating " << this->ix << "x" << this->iy << "px image...\t\t";

    this->planegrid_log =  new float[this->ix * this->iy];
    this->planegrid_real =  new float[this->ix * this->iy];

    // create a canvas
    this->plt = new Plotter(this->ix, this->iy);
    for(int i=0; i<this->ix; i++) {
        for(int j=0; j<this->iy; j++) {
            float x = _v1[0] * float(i - this->ix / 2) / _scale + _v2[0] * float(j - this->iy / 2) / _scale + _s[0];
            float y = _v1[1] * float(i - this->ix / 2) / _scale + _v2[1] * float(j - this->iy / 2) / _scale + _s[1];
            float z = _v1[2] * float(i - this->ix / 2) / _scale + _v2[2] * float(j - this->iy / 2) / _scale + _s[2];
            float val = this->sf->get_value_interp(x,y,z);
            if(negative_values) {
                if(val < -10) {
                    this->planegrid_log[j * this->ix + i] = -log10(-val);
                } else if(val > 10) {
                    this->planegrid_log[j * this->ix + i] = log10(val);
                } else {
                    this->planegrid_log[j * this->ix + i] = val / 10.0;
                }
            } else {
                this->planegrid_log[j * this->ix + i] = log10(val);
            }
            this->planegrid_real[j * this->ix + i] = val;
        }
    }

    std::cout << "[Done]" << std::endl;
}

void PlaneProjector::isolines(unsigned int bins) {
    float binsize = (this->max - this->min) / float(bins + 1);
    for(float val = this->min; val < this->max; val += binsize) {
        if(val < -1) {
            this->draw_isoline(-pow(10,-val));
        }
        if(val > 1) {
            this->draw_isoline(pow(10,val));
        }
    }
    this->draw_isoline(0);
}

void PlaneProjector::draw_isoline(float val) {
    for(unsigned int i=1; i<uint(this->ix-1); i++) {
        for(unsigned int j=1; j<uint(this->iy-1); j++) {
            if(this->is_crossing(i,j,val)) {
                this->plt->draw_filled_rectangle(i,j, 1, 1, Color(0,0,0));
            }
        }
    }
}

bool PlaneProjector::is_crossing(const unsigned int &i, const unsigned int &j, const float &val) {
    if(this->planegrid_real[(j-1) * this->ix + (i)] < val && this->planegrid_real[(j+1) * this->ix + (i)] > val) {
        return true;
    }
    if(this->planegrid_real[(j-1) * this->ix + (i)] > val && this->planegrid_real[(j+1) * this->ix + (i)] < val) {
        return true;
    }
    if(this->planegrid_real[(j) * this->ix + (i-1)] > val && this->planegrid_real[(j) * this->ix + (i+1)] < val) {
        return true;
    }
    if(this->planegrid_real[(j) * this->ix + (i-1)] < val && this->planegrid_real[(j) * this->ix + (i+1)] > val) {
        return true;
    }
    return false;
}

void PlaneProjector::plot() {
    for(unsigned int i=0; i<uint(this->ix); i++) {
        for(unsigned int j=0; j<uint(this->iy); j++) {
            this->plt->draw_filled_rectangle(i,j, 1, 1,
                this->scheme->get_color(this->planegrid_log[j * this->ix + i]));
        }
    }
}

void PlaneProjector::write(std::string filename) {
    plt->write(filename.c_str());
}

PlaneProjector::~PlaneProjector() {
    delete[] this->planegrid_log;
    delete[] this->planegrid_real;
}
