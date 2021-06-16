/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vegetationPBM

Description
	uncoupled solver for biomass particle thermo-chemical degradation
	
\*---------------------------------------------------------------------------*/
#include <iostream>
#include "argList.H"
#include "IOmanip.H"

using namespace std;
using namespace Foam;

#include "./particle_1D/particle_1D.h"


//A function that writes particle data to a file
void writeLine(FILE* outfile, double time, double numCells, std::vector<double> data)
{
    fprintf(outfile, "%.6f", time);
    for (int i = 0; i < numCells; i++)
        fprintf(outfile, "\t%.6f", data[i]);
    fprintf(outfile, "\n");
}


int main(int argc, char *argv[])
{

//inputs from gas phase
scalar localAirTemp = 300;           //current local gas temperature [K]
scalar localAirVel  = 0.3;           //current local gas velocity [m/s]
scalar Irradation   = 50.0e3;        //current incident radiation flux [W/m2]
scalar localO2Frac  = 0.0;           //current local oxygen vol fraction [-]

//particle geometry
word    geometry    = "rectangle";  //rectangle-cylinder-sphere
scalar  delta_i     = 1.25e-2;      //initial radius or half thickness [m]
scalar  length      = 1.0;          //rectangle/cylinder length [m]
scalar  area        = 1.0;          //rectangle area [m2]
scalar  resolution  = 100e-6;       //particle mesh resolution [m]

//particle material properties
scalar rho_m  = 1000;       // Mass density of liquid water [kg/m3]
scalar rho_vs = 663;        // Mass density of virgin solid [kg/m3]
scalar rho_c  = 132.6;      // Mass density of solid char [kg/m3]

scalar k_m    = 0.6;        // Conductivity of liquid water [W/m/K]
scalar k_vs   = 0.126;      // Conductivity of virgin solid [W/m/K]
scalar k_c    = 0.126;      // Conductivity of solid char [W/m/K]

scalar c_m    = 4.18e+3;    // Heat capacity of liquid water [J/kg/K]
scalar c_vs   = 2.52e+3;    // Heat capacity of virgin solid [J/kg/K]
scalar c_c    = 2.52e+3;   // Heat capacity of solid char [J/kg/K]

scalar eps_m  = 0.9;        // Surface emissivity of liquid water [-]
scalar eps_vs = 0.9;        // Surface emissivity of virgin solid [-]
scalar eps_c  = 0.9;        // Surface emissivity of solid char [-]

scalar eta_c  = 0.2;        // char yield [-]

//particle initial condition
scalar FMC    = 0.1;
scalar x_vs_i = 1.0/(1.0+(FMC*rho_vs/rho_m));  
scalar x_m_i  = 1.0-x_vs_i;                  
scalar T_i    = 300;

std::vector<double> T_list, x_m_list, x_vs_list;
T_list.assign(1,T_i);
x_m_list.assign(1,x_m_i);
x_vs_list.assign(1,x_vs_i);

//Simulation duration [s]
scalar simulationTime   = 1800;      
scalar globalTimeStep   = 1.0;
scalar globalTime       = 0.0;

//burning flag (0 inactive, 1 active)
bool status = true;

//output files headers
remove("particleTemporalOut.csv");
std::string header1 = "time (s), delta (m), mass (kg), vol (m3), Tsurf (k), MLR(kg/m3/s)";
ofstream writeTemporal("particleTemporalOut.csv", ios::app);
writeTemporal << header1 << endl;

remove("temperature.dat");
FILE* TempFile;
TempFile = fopen("temperature.dat", "w+");

remove("moisture.dat");
FILE* MoistureFile;
MoistureFile = fopen("moisture.dat", "w+");

remove("solid.dat");
FILE* SolidFile;
SolidFile = fopen("solid.dat", "w+");

remove("cellCenter.dat");
FILE* CoordFile;
CoordFile = fopen("cellCenter.dat", "w+");

//output fields
scalar  particleDelta, particleMass, particleVol;
scalar  T_surf, h_conv;
scalar  volProdRate_tot = 0.0;      // [kg/m3/s]
scalar  volProdRate_fuel = 0.0;     // [kg/m3/s]
scalar  volProdRate_H2O = 0.0;      // [kg/m3/s]
scalar  volHRR = 0.0;               // [J/m3/s]


// --------------------------------- run for a certain global time -----------------------------------

//initialize
particle_1D thisParticle;

thisParticle.set(
                    geometry,
                    resolution, 
                    delta_i, 
                    length, 
                    area,
                    eta_c,
                    rho_m,
                    rho_vs,
                    rho_c, 
                    c_m,
                    c_vs,
                    c_c, 
                    k_m,
                    k_vs,
                    k_c, 
                    eps_m,
                    eps_vs,
                    eps_c, 
                    T_list, 
                    x_m_list, 
                    x_vs_list
                );

double particleMass_i = thisParticle.getMass();

//time loop
while ((globalTime < simulationTime) && status)
{
    globalTime += globalTimeStep;

    thisParticle.stepForward(
                                globalTimeStep, 
                                localAirTemp, 
                                localAirVel, 
                                Irradation, 
                                localO2Frac,
                                status
                            );

    status = thisParticle.getState();
    particleDelta = thisParticle.getDelta();
    particleMass = thisParticle.getMass();
    particleVol = thisParticle.getVol();

    T_list      = thisParticle.getT();
    x_m_list    = thisParticle.getXm();
    x_vs_list   = thisParticle.getXvs();

    T_surf = thisParticle.getT().back();
    h_conv = thisParticle.getHconv();

    volProdRate_tot     = thisParticle.getMLR();
    volProdRate_fuel    = thisParticle.getGFRR(); 
    volProdRate_H2O     = volProdRate_tot - volProdRate_fuel;
    volHRR              = thisParticle.getCharHRR();

    writeTemporal << globalTime << "," << particleDelta << "," << particleMass << "," << particleVol 
    << "," << T_surf << "," << volProdRate_tot << endl;

    writeLine(TempFile, globalTime, thisParticle.getNcells(), T_list);
    writeLine(MoistureFile, globalTime, thisParticle.getNcells(), x_m_list);
    writeLine(SolidFile, globalTime, thisParticle.getNcells(), x_vs_list);
    writeLine(CoordFile, globalTime, thisParticle.getNcells(), thisParticle.getCellCenters());

    cout << "Global time (s): " << globalTime << ",  particleDelta: " << particleDelta << endl;

}


// ----------------------------------- check mass released --------------------------------------

// Initial mass of moisture and virgin solid (t = 0) [kg]
scalar Mm_i=0, Mvs_i=0, Mg_f=0;
const scalar pi = 3.14159265359;

if (geometry=="rectangle")
{
    Mm_i  = FMC*rho_vs*2*delta_i*area/(1+(FMC*rho_vs/rho_m));   // Mass of m
    Mvs_i =     rho_vs*2*delta_i*area/(1+(FMC*rho_vs/rho_m));   // Mass of vs 
}
else if (geometry=="cylinder")
{
    Mm_i  = FMC*rho_vs*pi*Foam::pow(delta_i,2.0)*length/(1+(FMC*rho_vs/rho_m)); // Mass of m
    Mvs_i =     rho_vs*pi*Foam::pow(delta_i,2.0)*length/(1+(FMC*rho_vs/rho_m)); // Mass of vs 
}

else if (geometry=="sphere")
{
    Mm_i  = FMC*rho_vs*4/3*pi*Foam::pow(delta_i,3.0)/(1+(FMC*rho_vs/rho_m)); // Mass of m
    Mvs_i =     rho_vs*4/3*pi*Foam::pow(delta_i,3.0)/(1+(FMC*rho_vs/rho_m)); // Mass of vs  
}


// Total mass of water vapor and volatiles released (t = inf) [kg]
if(eta_c == 0)      // Non-charring material
{
    Mg_f = Mm_i + Mvs_i;
}
else if(localO2Frac != 0)   // Charring material with char oxidation
{
    Mg_f = Mm_i + Mvs_i;
}
else                // Charring material with no char oxidation
{
    Mg_f = Mm_i + (1-eta_c)*Mvs_i;
}


cout << "Total mass of vapor/volatiles released (theory) = " << Mg_f << endl;
cout << "Total mass of vapor/volatiles released (simulation) = " << particleMass_i - particleMass << endl;


return 0;

}





