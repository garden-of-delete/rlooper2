//
// Created by Robert Stolz on 6/27/17.
//
#ifndef RLOOPER2_BIOPHYSICS_H
#define RLOOPER2_BIOPHYSICS_H

#include <math.h>

//Define global constants
//#define k_b 1.3806503e-23
#define R 0.0019858775

#define pi 3.14159265359

//globally relevant biophysics related functions.
long double compute_boltzmann_factor(double E, double T){
    return exp(-1*E/(R*T));
}


#endif //RLOOPER2_BIOPHYSICS_H
