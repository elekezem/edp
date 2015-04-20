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


PlaneProjector::PlaneProjector(ScalarField* _sf) {
    this->scheme = new ColorScheme(-5,5);
    this->sf = _sf;
}

void PlaneProjector::plot(Vector _v1, Vector _v2, Vector _s, float _scale, float li, float hi, float lj, float hj) {

    float _ix = int((hi - li) * _scale);
    float _iy = int((hj - lj) * _scale);

    std::cout << "Creating " << _ix << "x" << _iy << "px image." << std::endl;

    // create a canvas
    Plotter plt(_ix, _iy);
    for(float i=0; i<_ix; i++) {
        for(float j=0; j<_iy; j++) {
          float x = _v1[0] * float(i - _ix / 2) / _scale + _v2[0] * float(j - _iy / 2) / _scale + _s[0];
          float y = _v1[1] * float(i - _ix / 2) / _scale + _v2[1] * float(j - _iy / 2) / _scale + _s[1];
          float z = _v1[2] * float(i - _ix / 2) / _scale + _v2[2] * float(j - _iy / 2) / _scale + _s[2];
          float val = this->sf->get_value_interp(x,y,z);
          // check if there is a real value
          if(val == -1) {
            plt.draw_filled_rectangle(i,j, 1, 1, Color(0,0,0));
          } else {
            plt.draw_filled_rectangle(i,j, 1, 1, this->scheme->get_color(val));
          }
        }
    }
    plt.write("test.png");
}
