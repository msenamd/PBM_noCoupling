#ifndef PARTICLE_1D_H
#define PARTICLE_1D_H


using namespace std;
#include <vector>
#include <string>
#include <numeric>      // std::accumulate

#include "VegReaction.h"
#include "calcConvection.h"
#include "TDMA.h"

class particle_1D
{
public:
        particle_1D();
        virtual ~particle_1D();

    //member functions

   		//sets thermophysical properties of the particle
    	void set(
    				std::string geometry,
    				double resolution, 
    				double radius, 
    				double length, 
    				double area, 
    				double charYield,
    				double moistureDensity,
    				double solidDensity,
    				double charDensity, 
    				double moistureCapacity,
    				double solidCapacity,
    				double charCapacity, 
    				double moistureConductivity,
    				double solidConductivity,
    				double charConductivity, 
    				double moistureEmissivity,
    				double solidEmissivity,
    				double charEmissivity, 
					std::vector<double> currentTemp, 
					std::vector<double> current_x_m, 
					std::vector<double> current_x_vs
				);

    	//calculate minimum chemical timescale from used reaction mechanism
    	double tauChemical(); 

    	//calculate surface to volume ratio
    	double getSurfaceAreaToVolumeRatio();

		//clean up particle and destroy any dynamically allocated memory
		void destroy();

    	//step forward one (global) time step give local gas temperature, velocity, radiant flux and oxygen fraction
		void stepForward(const double globalTimeStep, const double T_g, const double u_g,
							const double G , const double x_O2_g, bool active);

protected:
    const double pi = 3.14159265358979323846;
	const double sigma = 5.67e-8;     //% Stefan-Boltzmann constant [W/m2/K4]
	const double R = 8.3144598;           //Gas constant (J/(mol*K))

private:
		
    	std::string shape;      			// particle geometry
    	double delta_i;                		// Initial half-thickness (rectangle) or radius (cyl,sphere) of the particle
    	double areaRectangle;            	// Area of exposed surface (for rectangular particles only) (m2)
    	double lengthCylinder;             	// Cylinder length (for cylindrical particles only) (m)
		
		double eta_c;                  		// Char yield (-)

		double rho_m;                  		// Moisture Density (kg/m3)
		double rho_vs;                 		// Solid Density (kg/m3)
		double rho_c;                  		// Char Density (kg/m3)

		double c_m;                  		// Moisture Capacity (J/kg K)
		double c_vs;                 		// Solid Capacity (J/kg K)
		double c_c;                  		// Char Capacity (J/kg K)

		double k_m;                  		// Moisture Conductivity (W/mk)
		double k_vs;                 		// Solid Conductivity (W/mk)
		double k_c;                  		// Char Conductivity (W/mk)	

		double e_m;                  		// Moisture Emissivity (-)
		double e_vs;                 		// Solid Emissivity (-)
		double e_c;                  		// Char Emissivity (-)	

		double delta;                  		//instantaneous particle size (radius or half-thickness)
    	double Tsurf;                  		//temperature of the surface (K)
    	double Tcore;                  		//temperature of the inner most cell (K)

		double MLR;                    		//particle volumetric mass loss rate (including moisture and fuel) (kg/m3.s)
		double GFRR;				   		//gaseous fuel volumetric production rate (kg/m3.s)
		double charHRR;                		//char oxid. volumetric heat release rate (J/m3.s)

		bool state;

    	//reaction mechanism and additional species
    	VegReaction R1;
    	VegReaction R2;
    	VegReaction R3;

		int nx_i;                       	// initial number of grid cells (-)
		int nx;                         	// number of cells (-)
		int m;                          	// local time iteration counter (s)
		int m_f;                        	// total number of particle local time-steps (-)
		double dx_i;                    	// Initial grid cell size (-)
        double dt;                      	// integration step size (s)
		double h_conv;                  	// convective heat transfer coefficient (W/m2K)
		double volume_i;                 	// particle initial volume (m3)
		double volume;						// particle final volume (m3)
		double eps_surf;                	// initial particle surface emissivity (-)
		double q_surf;                  	// Net surface heat flux (W/m2)
		double kp_surf;                 	// Surface thermal conductivity (W/mK)
		double dx_surf;                 	// grid spacing at the surface (m)
		double localTime;               	// local particle time, tarts from 0 (s)
		double mass_p;						// particle mass (kg)

        double FO_L;                    	// intermediate variable used in energy equation discretization
        double FO_R;                    	// intermediate variable used in energy equation discretization

		std::vector<double> Temp;            // particle temperature (K)
		std::vector<double> x_m, x_vs, x_c;  // volume fractions of moisture, virgin solid and char (-)
		std::vector<double> dx;              // cell size (length of rectlinear grid or radius of spherical and cylindrical shells) (m)
		std::vector<double> dV;              // cell volume (m3)
		std::vector<double> rho_dV_cell;     // cell mass (kg)
		std::vector<double> xRight;          // location of downwind face (m)
		std::vector<double> xCenter;         // location of cell center (m)
		std::vector<double> MLR_cell;        // Mass Loss Rate of each cell (kg/s)
		std::vector<double> GFRR_cell;       // Gaseous fuel release rate from each cell (kg/s)
		std::vector<double> charHRR_cell;    // Char oxid. heat release rate from each cell (kg/s)

        std::vector<double> S_pos;           //downwind face area (m2)
        std::vector<double> S_neg;           //upwind face area (m2)
        std::vector<double> kp;              //particle thermal conductivity (W/mK)
        std::vector<double> rho_times_cp;	 //rho*cp
        std::vector<double> Qdotp;           //volumetric rate of heat production/consumption (W/m3)

		//variable declaration at old local timestep
		int nx_old;
		double dt_old;
		std::vector<double> xCenter_old;
		std::vector<double> xRight_old;
		std::vector<double> dV_old;
		std::vector<double> dx_old;
		std::vector<double> Temp_old;
		std::vector<double> x_m_old;
		std::vector<double> x_vs_old;
		std::vector<double> x_c_old;
		double q_surf_old;
		double Tsurf_old;

		//dummy std::vectors for TDMA conversion
		std::vector<double> vectA;
		std::vector<double> vectB;
		std::vector<double> vectC;
		std::vector<double> vectD;

		// operator splitting method sub-step particle temperature
		std::vector<double> Temp0, Temp1, Temp2, Temp3; 

        //calculate time step size based on 0.5 Fourier number
        double tauDiff(std::vector<double> currentTemp, std::vector<double> current_x_m, std::vector<double> current_x_vs, double deltaX);

		// Set computational grid (works for rectlinear, cylindrical or spherical grids)
		tuple<std::vector<double>, std::vector<double>, std::vector<double>> setMesh(
			int& numCells, std::vector<double>& sizeCell,
			const std::string& gridGeom, const double& Area, const double& Length);

		/* reconstructs the grid (in case of volume change)
		   - sets the first three std::vectors as: downwind face coord - cell center coord - cell spacing */
		void moveMesh(std::vector<double>& xR, std::vector<double>& xC, std::vector<double>& deltaX,
			const int& numCells, std::vector<double>& vol, const std::string& gridGeom , const double& Area, const double& Length);


		//solve mass conservation
		void mass_conservation(double x_O2_g);

		//solve energy conservation - explicit
		std::vector<double> energy_conservation_explicit(double x_O2_g);

		//solve energy conservation - reaction step
		std::vector<double> energy_conservation_reactionstep(double x_O2_g);

		//solve energy conservation - diffusion step
		void energy_conservation_diffusionstep(double T_g, double G, double x_O2_g);

		//calculate total particle volume
		double getVolume(double);

        //intermediate calculation for energy equation, for left face
        void set_FO_L(const int &i , const double &dt_);

        //intermediate calculation for energy equation, for right face
        void set_FO_R(const int &i , const double &dt_);

		//adjust the local time-step of particle integration
		void adjustTimeStep(double globalTimeStep, double x_O2_g);

		//update particle surface temperature
		void updateSurface(double T_g, double G);

		/* Remeshing
			-Reduce nuber of cells while particle is burning
			-Maintain uniform grid with spacing close to initial value */
		void remeshing();
		void interpolateOnNewMesh(int nx_new);
		std::vector<double> solveCostFunction(std::vector<double> q0, std::vector<double> q_old, int nx_new, int nx_old);

		//cap vector to min and max values
		void checkBounds(std::vector<double> q);

public:
		//Access
		bool getState();							//return particle stat (-)
		double getNcells();							//return number of cells (-)
		double getDelta();							//return particle radius or half thickness (m)
		double getMass();							//return particle mass (kg)
		double getVol();							//return particle volume (m3)
		double getMLR();							//return volumetric mass lss rate (kg/m3/s)
		double getGFRR();							//return volumetric production of gaseous mixtures (kg/m3/s)
		double getCharHRR();						//return volumetric heat release rate (J/m3/s)
		double getHconv();							//return heat transfer coeff. (W/m2K)
		std::vector<double> getCellCenters();		//return list of location of cell centers (m)
		std::vector<double> getT();					//return list of cell center temperatures (K)
		std::vector<double> getXm();				//return list of cell center moisture vol. frac. (-)
		std::vector<double> getXvs();				//return list of cell center virgin solid vol. frac. (-)

};

#endif // PARTICLE_1D_H
