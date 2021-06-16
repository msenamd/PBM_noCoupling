/*
 * calcConvection.h
 * This function calculates heat transfer coefficient from empirical
 * correlations of rectangle, cylinder & sphere
 * Inputs: gas temperature, surface temp, gas velocity, particle size, particle geometry
 *
 * rectangle and cylinder: Ref: Churchill and Bernstein
 * sphere: Ref: Whitaker
 *
 *  Created on: October 9, 2019
 *      Author: mmahmed
 */
#include <iostream>
#include <vector>
#include <cmath>
#include "particle_1D.h"

using namespace std;

#ifndef CALCCONVECTION_H
#define CALCCONVECTION_H

double get_h(double,double,double,double,std::string);

#endif // CALCCONVECTION_H
