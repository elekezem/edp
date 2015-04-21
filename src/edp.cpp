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
#include <tclap/CmdLine.h> // parsing command line arguments
#include "mathtools.h"
#include "scalar_field.h"
#include "planeprojector.h"

int main(int argc, char *argv[]) {
    // command line grabbing
    try {
        TCLAP::CmdLine cmd("Projects the electrondensity of a CHGCAR file onto a image file.", ' ', "0.9");

        //**************************************
        // declare values to be parsed
        //**************************************
        TCLAP::ValueArg<std::string> arg_name("o","filename","Filename to print to",true,"test.png","string");
        cmd.add(arg_name);
        TCLAP::ValueArg<std::string> arg_sp("p","starting_point","Start point of cutting plane",true,"(0.5,0.5,0.5)","3d-vector");
        cmd.add(arg_sp);
        TCLAP::ValueArg<std::string> arg_v("v","vector1","Plane Vector 1",true,"(1,0,0)","3d-vector");
        cmd.add(arg_v);
        TCLAP::ValueArg<std::string> arg_w("w","vector2","Plane Vector 2",true,"(0,0,1)","3d-vector");
        cmd.add(arg_w);
        TCLAP::ValueArg<unsigned int> arg_s("s","scale","Scaling in px/angstrom",true, 200,"unsigned integer");
        cmd.add(arg_s);
        TCLAP::SwitchArg arg_negative("n","negative_values","CHGCAR can contain negative values", cmd, false);

        cmd.parse(argc, argv);

        //**************************************
        // parsing values
        //**************************************
        std::string filename = arg_name.getValue();

        pcrecpp::RE re("^([0-9.-]+),([0-9.-]+),([0-9.-]+)$");
        std::string sp = arg_sp.getValue();
        float sp_in[3];
        re.FullMatch(sp.c_str() , &sp_in[0], &sp_in[1], &sp_in[2]);
        Vector s(sp_in[0],sp_in[1],sp_in[2]);

        std::string v = arg_v.getValue();
        float v_in[3];
        re.FullMatch(v.c_str() , &v_in[0], &v_in[1], &v_in[2]);
        Vector v1(v_in[0],v_in[1],v_in[2]);

        std::string w = arg_w.getValue();
        float w_in[3];
        re.FullMatch(w.c_str() , &w_in[0], &w_in[1], &w_in[2]);
        Vector v2(w_in[0],w_in[1],w_in[2]);

        float scale = arg_s.getValue();

        bool negative_values = arg_negative.getValue();

        //**************************************
        // start running the program
        //**************************************
        std::cout << "Running EDP" << std::endl;

        // read in field
        ScalarField sf("CHGCAR");
        sf.read(true);

        // define intervals in Angstrom
        float interval = 20.0;
        float li = -interval;
        float hi = interval;
        float lj = -interval;
        float hj = interval;

        float color_interval = 5;

        PlaneProjector pp(&sf, -color_interval, color_interval);
        pp.extract(v1, v2, s, scale, li, hi, lj, hj, negative_values);
        pp.plot();
        pp.isolines(int(color_interval + 1)*2, negative_values);
        pp.write(filename);

        return 0;
    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
