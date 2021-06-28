#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <vector>
#include <ctime>
#include <cstring>
#include <tuple>
#include <omp.h>

#include "./particle_1D/particle_1D.h"

using namespace std;

//A function that writes particle data to a file
void writeLine(FILE* outfile, double time, double numCells, std::vector<double> data)
{
    fprintf(outfile, "%.6f", time);
    for (int i = 0; i < numCells; i++)
        fprintf(outfile, "\t%.6f", data[i]);
    fprintf(outfile, "\n");
}


int main()
{

//inputs from gas phase
double localAirTemp = 300;           //current local gas temperature [K]
double localAirVel  = 0.3;           //current local gas velocity [m/s]
double Irradation   = 50.0e3;        //current incident radiation flux [W/m2]
double localO2Frac  = 0.0;           //current local oxygen vol fraction [-]

//particle geometry
string    geometry    = "rectangle";  //rectangle-cylinder-sphere
double  delta_i     = 1.25e-2;      //initial radius or half thickness [m]
double  length      = 1.0;          //rectangle/cylinder length [m]
double  area        = 1.0;          //rectangle area [m2]
double  resolution  = 100e-6;       //particle mesh resolution [m]

//particle material properties
double rho_m  = 1000;       // Mass density of liquid water [kg/m3]
double rho_vs = 663;        // Mass density of virgin solid [kg/m3]
double rho_c  = 132.6;      // Mass density of solid char [kg/m3]

double k_m    = 0.6;        // Conductivity of liquid water [W/m/K]
double k_vs   = 0.126;      // Conductivity of virgin solid [W/m/K]
double k_c    = 0.126;      // Conductivity of solid char [W/m/K]

double c_m    = 4.18e+3;    // Heat capacity of liquid water [J/kg/K]
double c_vs   = 2.52e+3;    // Heat capacity of virgin solid [J/kg/K]
double c_c    = 2.52e+3;   // Heat capacity of solid char [J/kg/K]

double eps_m  = 0.9;        // Surface emissivity of liquid water [-]
double eps_vs = 0.9;        // Surface emissivity of virgin solid [-]
double eps_c  = 0.9;        // Surface emissivity of solid char [-]

double eta_c  = 0.2;        // char yield [-]

//particle initial condition
double FMC    = 0.1;
double x_vs_i = 1.0/(1.0+(FMC*rho_vs/rho_m));  
double x_m_i  = 1.0-x_vs_i;                  
double T_i    = 300;

std::vector<double> T_list, x_m_list, x_vs_list;
T_list.assign(1,T_i);
x_m_list.assign(1,x_m_i);
x_vs_list.assign(1,x_vs_i);

//Simulation duration [s]
double simulationTime   = 100;      
double globalTimeStep   = 1.0;
double globalTime       = 0.0;

//burning flag (0 inactive, 1 active)
bool status = true;

//output files headers
remove("particleTemporalOut.csv");
std::string header1 = "time (s), delta (m), mass (kg), vol (m3), Tsurf (k), Tcore (k), MLR(kg/m3/s)";
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
double  particleDelta, particleMass, particleVol;
double  T_surf, T_core, h_conv;
double  volProdRate_tot = 0.0;      // [kg/m3/s]
double  volProdRate_fuel = 0.0;     // [kg/m3/s]
double  volProdRate_H2O = 0.0;      // [kg/m3/s]
double  volHRR = 0.0;               // [J/m3/s]


// --------------------------------- run for a certain global time -----------------------------------

std::clock_t c_start = std::clock();

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
    T_core = thisParticle.getT().front();

    h_conv = thisParticle.getHconv();

    volProdRate_tot     = thisParticle.getMLR();
    volProdRate_fuel    = thisParticle.getGFRR(); 
    volProdRate_H2O     = volProdRate_tot - volProdRate_fuel;
    volHRR              = thisParticle.getCharHRR();

    writeTemporal << globalTime << "," << particleDelta << "," << particleMass << "," << particleVol 
    << "," << T_surf << "," << T_core << "," << volProdRate_tot << endl;

    writeLine(TempFile, globalTime, thisParticle.getNcells(), T_list);
    writeLine(MoistureFile, globalTime, thisParticle.getNcells(), x_m_list);
    writeLine(SolidFile, globalTime, thisParticle.getNcells(), x_vs_list);
    writeLine(CoordFile, globalTime, thisParticle.getNcells(), thisParticle.getCellCenters());

    cout << "Global time (s): " << globalTime << ",  particleDelta: " << particleDelta << endl;

}


// ----------------------------------- check mass released --------------------------------------

// Initial mass of moisture and virgin solid (t = 0) [kg]
double Mm_i=0, Mvs_i=0, Mg_f=0;
const double pi = 3.14159265359;

if (geometry=="rectangle")
{
    Mm_i  = FMC*rho_vs*2*delta_i*area/(1+(FMC*rho_vs/rho_m));   // Mass of m
    Mvs_i =     rho_vs*2*delta_i*area/(1+(FMC*rho_vs/rho_m));   // Mass of vs 
}
else if (geometry=="cylinder")
{
    Mm_i  = FMC*rho_vs*pi*std::pow(delta_i,2.0)*length/(1+(FMC*rho_vs/rho_m)); // Mass of m
    Mvs_i =     rho_vs*pi*std::pow(delta_i,2.0)*length/(1+(FMC*rho_vs/rho_m)); // Mass of vs 
}

else if (geometry=="sphere")
{
    Mm_i  = FMC*rho_vs*4/3*pi*std::pow(delta_i,3.0)/(1+(FMC*rho_vs/rho_m)); // Mass of m
    Mvs_i =     rho_vs*4/3*pi*std::pow(delta_i,3.0)/(1+(FMC*rho_vs/rho_m)); // Mass of vs  
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


std::clock_t c_end = std::clock();

double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

return 0;

}





