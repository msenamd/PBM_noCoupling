/*
 * TDMA.h
 * This is a tri-diagonal matrix algorithm (Thomas Algorithm)
 * arrays a,b,c,d are inputs (not changed by the code)
 * the output array f is the solution
 *
 *  Created on: October 14, 2019
 *      Author: mmahmed
 */

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#ifndef TDMA
#define TDMA
 std::vector<double> solveTriDiag(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
#endif
