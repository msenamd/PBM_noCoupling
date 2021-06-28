#include "particle_1D.h"
#include "energy_conservation_explicit.h"
#include "energy_conservation_reactionstep.h"
#include "energy_conservation_diffusionstep.h"
#include "mass_conservation.h"
#include "remeshing.h"
#include "particle_surface.h"

particle_1D::particle_1D()
{
	shape			= " ";
    delta_i  		= 0.0;
    areaRectangle 	= 0.0;
    lengthCylinder  = 0.0;

	eta_c 		  	= 0.0;

	rho_m			= 0.0;
	rho_vs			= 0.0;
	rho_c			= 0.0;

	c_m				= 0.0;
	c_vs			= 0.0;
	c_c				= 0.0;

	k_m				= 0.0;
	k_vs			= 0.0;
	k_c				= 0.0;

	e_m				= 0.0;
	e_vs			= 0.0;
	e_c				= 0.0;

	delta 			= 0.0;
	Tsurf 			= 0.0;
	Tcore 			= 0.0;

	MLR 			= 0.0;
	GFRR 			= 0.0;
	charHRR 		= 0.0;

	state 			= true;

    R1.set_knownReaction("drying");
    R2.set_knownReaction("pyrolysis");
    R3.set_knownReaction("char oxidation");

	nx_i 			= 0;
	nx 				= 0;
	m_f				= 0;
	m 				= 0;
	dx_i 			= 0.0;
	dt 				= 0.0;
	h_conv 			= 0.0;
	volume_i 		= 0.0;
	volume 			= 0.0;
	eps_surf 		= 0.0;
	q_surf 			= 0.0;
	kp_surf 		= 0.0;
	dx_surf 		= 0.0;
	localTime		= 0.0;
	mass_p			= 0.0;

	FO_L 			= 0.0;
	FO_R 			= 0.0;

	nx_old 			= 0;
	dt_old 			= 0.0;
	q_surf_old 		= 0.0;
	Tsurf_old 		= 0.0;
}

particle_1D::~particle_1D()
{
	destroy();
}


// Calculate minimum chemical time-scale
double particle_1D::tauChemical(){

    double tau_R1_min,tau_R2_min;

    tau_R1_min=1/(R1.get_A()*exp(-R1.get_Ta()/500)); //500 is maximum temperature range of R1
    tau_R2_min=1/(R2.get_A()*exp(-R2.get_Ta()/1000));//1000 is maximum temperature range of R2

    return min(tau_R1_min,tau_R2_min);
}


// Sets geometry, dimensions and properties
void particle_1D::set
(
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
)
{
    shape 	= geometry;
    delta_i = radius;
    lengthCylinder 	= length;
    areaRectangle 	= area;

    eta_c 	= charYield;

    rho_m  	= moistureDensity;
    rho_vs 	= solidDensity;
    rho_c  	= charDensity;

    c_m  	= moistureCapacity;
    c_vs 	= solidCapacity;
    c_c  	= charCapacity;

    k_m  	= moistureConductivity;
    k_vs 	= solidConductivity;
    k_c  	= charConductivity;

    e_m  	= moistureEmissivity;
    e_vs 	= solidEmissivity;
    e_c  	= charEmissivity;

	nx_i 	= round(delta_i / (resolution)); 	//spatial resolution ~ 10-microns
	nx_i 	= max(nx_i, 5);          			//initial number of cells, keep min of 5 cells

    if(currentTemp.size() == 1)
	{
		nx = nx_i;
		dx_i = delta_i / nx;
		dx.assign(nx, dx_i);
		Temp.assign(nx, currentTemp.front());
		x_m.assign(nx, current_x_m.front());
		x_vs.assign(nx, current_x_vs.front());
		x_c.assign(nx, max(0.0, min(1.0, 1.0 - current_x_m.front() - current_x_vs.back())));
	}
	else
	{
		nx=currentTemp.size();
		dx_i = delta_i / nx;
		dx.assign(nx, dx_i);
		Temp = currentTemp;
		x_m = current_x_m;
		x_vs = current_x_vs;
		x_c.assign(x_m.size(),0.0);
	
		for(unsigned int j=0; j<x_m.size(); j++)
		{
			x_c[j] = 1.0-x_m[j]-x_vs[j];
		}
	}
	
	//set time step
	dt = min(tauChemical()*0.1,tauDiff(Temp,x_m,x_vs,dx_i));
	
	Tsurf = Temp.front();
	Tcore = Temp.back();
	
	xRight = get<0>(setMesh(nx, dx, shape, areaRectangle, lengthCylinder));
	xCenter = get<1>(setMesh(nx, dx, shape, areaRectangle, lengthCylinder));
	dV = get<2>(setMesh(nx, dx, shape, areaRectangle, lengthCylinder));

	delta = delta_i;
	volume_i = getVolume(delta_i);
	volume = volume_i;

	mass_p = 0.0;	
	for (int i = 0 ; i < nx_i ; i++)
	{
		mass_p += (rho_m*x_m[i] + rho_vs*x_vs[i] + rho_c*x_c[i]) * dV[i];
	}
	if(shape=="rectangle")
	{
		mass_p = 2*mass_p;
	}

	eps_surf = e_m * current_x_m.front()
				+ e_vs * current_x_vs.front()
				+ e_c * (1.0-current_x_m.front()-current_x_vs.front());		
}


/* Set computational grid (works for rectlinear, cylindrical or spherical grids)
   - sets the first three std::vectors as: downwind face coord - cell center coord - cell volume
   - overwrites the first three std::vectors to minimize memory usage */
tuple<std::vector<double>, std::vector<double>, std::vector<double>> particle_1D::setMesh(
	int& numCells, std::vector<double>& sizeCell,
	const std::string& gridGeom, const double& Area, const double& Length)
{
	std::vector<double> xR, xC, vol;
	xR.assign(numCells, 0.0);    //downwind face location
	xC.assign(numCells, 0.0);    //cell center location
	vol.assign(numCells, 0.0);   //cell volume
	double xL1 = 0.0;            //first upwind face location (always at 0.0)

	xR[0] = sizeCell[0];
	xC[0] = 0.5 * (xL1 + xR[0]);
	for (int i = 1; i < numCells; i++)
	{
		xR[i] = sizeCell[i] + xR[i-1];
		xC[i] = 0.5 * (xR[i-1] + xR[i]);

	}

	if (gridGeom == "rectangle")
	{
		for (int i = 0; i < numCells; i++)
		{
			vol[i] = sizeCell[i] * Area;
		}
	}
	else if (gridGeom == "cylinder")
	{
		vol[0] = pi * sizeCell[0] * (xR[0]) * Length;
		for (int i=1; i < numCells; i++)
		{
			vol[i] = pi * sizeCell[i] * (xR[i] + xR[i-1]) * Length;
		}
	}
	else if (gridGeom == "sphere")
	{
		vol[0] = 4.0/3.0 * pi * sizeCell[0] * pow(xR[0] , 2.0);
		for (int i = 1; i < numCells; i++)
		{
			vol[i] = 4.0/3.0 * pi * sizeCell[i] *
				( pow(xR[i] , 2.0) + xR[i] * xR[i-1] + pow(xR[i-1] , 2.0));
		}
	}
	return make_tuple(xR,xC,vol);
}

// Calculate time step size based on 0.5 Fourier number
double particle_1D::tauDiff(std::vector<double> currentTemp, std::vector<double> current_x_m, std::vector<double> current_x_vs, double deltaX)
{
    double dt_diff = 0.0;
    double FO = 0.5;      //Fourier Number
    double alpha_max = 0.0;
    double FMC_i;
    double x_vs_tmp, FMC_tmp, x_m_tmp, x_c_tmp, rhoc_p_tmp, k_p_tmp;
	
	double x_m_,x_vs_;
	
	x_m_=current_x_m.front();
	x_vs_=current_x_vs.front();

    FMC_i=(rho_m*x_m_)/(rho_vs*(x_vs_));

    for(int i=1;i<101;i++){
        x_vs_tmp = (i-1)/100;   // Cover range 0 <= x_vs_tmp <= 1
        for(int j=1;j<21;j++){
            FMC_tmp = (j-1)*FMC_i/20;   // Cover range 0 <= FMC_tmp <= FMC
            x_m_tmp = (rho_vs/rho_m)*FMC_tmp*x_vs_tmp;
            x_c_tmp = 1-x_vs_tmp-x_m_tmp;

            rhoc_p_tmp = rho_m*c_m*x_m_tmp
                        + rho_vs*c_vs*x_vs_tmp
                        + rho_c*c_c*x_c_tmp;

            k_p_tmp    = k_m*x_m_tmp + k_vs*x_vs_tmp
                        + k_c*x_c_tmp;

            alpha_max = max(alpha_max,(k_p_tmp/rhoc_p_tmp));
        }
    }

dt_diff =  FO*pow(deltaX,2.0)/alpha_max;
return dt_diff;
}


// Calculates surface area to volume rate
double particle_1D::getSurfaceAreaToVolumeRatio()
{
	double areaToVolume = 0.0;

	if (shape == "rectangle")
	{
		areaToVolume = 1.0 / (2.0 * xRight.back());
	}
	else if (shape == "cylinder")
	{
		areaToVolume = 2.0 / xRight.back();
	}
	else if (shape == "sphere")
	{
		areaToVolume = 3.0 / xRight.back();
	}
	return areaToVolume;
}

// Destroy any dynamically allocated memory
void particle_1D::destroy()
{

}

// Do one step in global time
void particle_1D::stepForward(const double globalTimeStep,
                                  const double T_g, const double u_g,
                                  const double G , const double x_O2_g, bool active)
{
	m = 0;
	dt = min(dt, globalTimeStep);  //in case of globalTimeStep>dt
	m_f = round(globalTimeStep / dt);

	if ((m_f * dt) != globalTimeStep) {    //in case of dt*n_f not equal globalTimeStep
		m_f = m_f + 1;
		dt = globalTimeStep / m_f;
	}

	h_conv = get_h(T_g, Tsurf, u_g, xRight.back(), shape);
	q_surf = eps_surf * G - eps_surf * sigma * pow(Tsurf, 4.0)
			+ h_conv * (T_g - Tsurf);

	localTime = 0.0;
	MLR = 0.0;
	GFRR = 0.0;
	charHRR = 0.0;

	state = active;

	// Local time loop
	while ((m < m_f) && (state == true))
	{
		m = m + 1;

		nx_old 			= nx;
		dV_old 			= dV;
		dx_old 			= dx;
		dt_old			= dt;
		Temp_old 		= Temp;
		x_m_old 		= x_m;
		x_vs_old 		= x_vs;
		x_c_old 		= x_c;
		q_surf_old 		= q_surf;
		Tsurf_old 		= Tsurf;

		adjustTimeStep(globalTimeStep, x_O2_g);

		localTime = localTime + dt;

	    // - Calculate temperature using operator splitting method (Strang's method, 3 steps)

	    	// Step 1: calculate temp1 as the temperature at sub-step t(n+0.25)
		    Temp0 = Temp_old;
		    Temp1 = energy_conservation_reactionstep(x_O2_g);

		    // Step 2: calculate temp2 as the temperature at sub-step t(n+0.75)
		    Temp0 = Temp1;
			vectA.assign(nx_old, 0.0);
			vectB.assign(nx_old, 0.0);
			vectC.assign(nx_old, 0.0);
			vectD.assign(nx_old, 0.0);
			
			energy_conservation_diffusionstep(T_g, G, x_O2_g);
			
			Temp2 = solveTriDiag(vectA,vectB,vectC,vectD);

			// Step 3: calculate temp3 as the temperature at sub-step t(n+1)
			Temp0 = Temp2;
			Temp3 = energy_conservation_reactionstep(x_O2_g);

			Temp = Temp3;

		// Solve mass conservation
		mass_conservation(x_O2_g);

		// Update particle size
		delta 	= xRight.back();
		volume 	= getVolume(delta);

		// Update computational grid (in case of volume change)
		moveMesh(xRight, xCenter, dx, nx_old, dV, shape, areaRectangle, lengthCylinder);

		// Calculate the net surface heat flux and the surface temperature
		kp_surf = k_m * x_m.back() + k_vs * x_vs.back() + k_c * x_c.back();
		eps_surf = e_m * x_m.back() + e_vs * x_vs.back() + e_c * x_c.back();
		dx_surf = xRight.back() - xCenter.back();
		h_conv = get_h(T_g, Tsurf_old, u_g, xRight.back(), shape);

		// Calculate the net surface heat flux and the surface temperature
		updateSurface(T_g, G);

		// Calculate volumetric mass loss rate and gas fuel released at this local time-step (for each cell)
		MLR_cell.assign(nx, 0.0);
		GFRR_cell.assign(nx, 0.0);
		charHRR_cell.assign(nx, 0.0);
		rho_dV_cell.assign(nx, 0.0);

		for (int i = 0; i < nx; i++)
		{
			rho_dV_cell[i] = dV[i] * (
							rho_m * x_m[i] 
							+ rho_vs * x_vs[i] 
							+ rho_c * x_c[i]
									 );

			MLR_cell[i] = dV[i] * (
				        rho_m * x_m[i] * R1.get_A() * exp(-R1.get_Ta() / Temp[i])
						+ rho_vs * x_vs[i] * R2.get_A() * exp(-R2.get_Ta() / Temp[i]) * (1 - eta_c)
						+ rho_c * x_c[i] * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp[i])
								 );

			GFRR_cell[i] = dV[i] * (
							rho_vs * x_vs[i] * R2.get_A() * exp(-R2.get_Ta() / Temp[i]) * (1 - eta_c)
							+ rho_c * x_c[i] * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp[i])
								 );

			charHRR_cell[i] = dV[i] * 
							rho_c * x_c[i] * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp[i]) * R3.get_deltaH();
		}

		// Update particle mass (kg); and volumetric outputs (int.dV/V)
		mass_p 	= accumulate(rho_dV_cell.begin(), rho_dV_cell.end(), 0.0f);	 	 			 //[kg]
		MLR 	= accumulate(MLR_cell.begin(), MLR_cell.end(), 0.0f) / volume;             //[kg/m3.s]
		GFRR 	= accumulate(GFRR_cell.begin(), GFRR_cell.end(), 0.0f) / volume;		     //[kg/m3.s]
		charHRR = accumulate(charHRR_cell.begin(), charHRR_cell.end(), 0.0f) / volume; 	 //[J/m3.s]

		//correction to account for full rectangualr particle
		if (shape == "rectangle")
		{
			mass_p = mass_p	 * 2;
			MLR = MLR * 2;
			GFRR = GFRR * 2;
			charHRR = charHRR * 2;
		}

		//update particle status
		if ((eta_c != 0.0) && (x_vs.front() < 0.001))
		{
			state = false;
			cout <<"state change to burned" << endl;
		}
		if ((eta_c == 0.0) && (delta/delta_i < 0.02))
		{
			state = false;
			cout <<"state change to burned" << endl;
		}
	
		//auto-remesh the particle if volume shrinks
		remeshing();

	} //end time loop
	
}


//intermediate calculation for energy equation, for left face
void particle_1D::set_FO_L(const int &i , const double &dt_)
{
    FO_L = 0.5 * (kp[i] + kp[i - 1]) / rho_times_cp[i] * S_neg[i] * dt_ /
        (dV[i] * (xCenter[i] - xCenter[i - 1]));
}

//intermediate calculation for energy equation, for right face
void particle_1D::set_FO_R(const int &i , const double &dt_)
{
    FO_R = 0.5 * (kp[i] + kp[i + 1]) / rho_times_cp[i] * S_pos[i] * dt_ /
        (dV[i] * (xCenter[i + 1] - xCenter[i]));
}

//reconstructs the grid(in case of volume change)
void particle_1D::moveMesh(std::vector<double>& xR, std::vector<double>& xC, std::vector<double>& deltaX,
	const int& numCells, std::vector<double>& vol, const std::string& gridGeom, const double& Area, const double& Length)
{
	if (gridGeom == "rectangle")
	{
		deltaX[0] = vol[0] / Area;
		xR[0] = deltaX[0];
		xC[0] = 0.5 * xR[0];
		for (int i = 1; i < numCells; i++)
		{
			deltaX[i] = vol[i] / Area;
			xR[i] = deltaX[i] + xR[i-1];
			xC[i] = 0.5 * (xR[i-1] + xR[i]);
		}
	}
	else if (gridGeom == "cylinder")
	{
		deltaX[0] = pow((vol[0] / Length / pi) , 0.5);
		xR[0] = deltaX[0];
		xC[0] = 0.5 * xR[0];
		for (int i = 1; i < numCells; i++)
		{
			deltaX[i] = pow(((vol[i] / Length + pi * pow(xR[i-1] , 2.0)) / pi) , 0.5) - xR[i-1];
			xR[i] = deltaX[i] + xR[i-1];
			xC[i] = 0.5 * (xR[i-1] + xR[i]);
		}
	}
	else if (gridGeom == "sphere")
	{
		deltaX[0] = pow((3.0 * vol[0] / 4.0 / pi), (1.0 / 3.0));
		xR[0] = deltaX[0];
		xC[0] = 0.5 * xR[0];
		for (int i = 1; i < numCells; i++)
		{
			deltaX[i] = pow( ((3.0*vol[i]/4.0/pi + pow(xR[i-1] , 3.0) )) , (1.0/3.0)) - xR[i - 1];
			xR[i] = deltaX[i] + xR[i - 1];
			xC[i] = 0.5 * (xR[i - 1] + xR[i]);
		}
	}
}

//calculate total particle volume
double particle_1D::getVolume(double XX)
{
	double Vp=0.0;
	if (shape == "rectangle") {
		Vp = 2 * XX * areaRectangle;
	}
	else if (shape == "cylinder")
	{
		Vp = pi * pow(XX, 2.0) * lengthCylinder;
	}
	else if (shape == "sphere")
	{
		Vp = 4.0 / 3.0 * pi * pow(XX, 3.0);
	}
	return Vp;
}

// Adjusting time-step
void particle_1D::adjustTimeStep(double globalTimeStep, double x_O2_g)
{
	// -Estimate Tp(n+1) from explicit treatment
	Temp = energy_conservation_explicit(x_O2_g);

	// -Calculate maximum temperature difference
	vector<double> DeltaTemp;
	DeltaTemp.assign(nx_old, 0.0);
	for (int i = 0; i < nx; i++)
	{
		DeltaTemp[i] = abs(Temp[i] - Temp_old[i]);
	}
	double DeltaTemp_max = *max_element(DeltaTemp.begin(), DeltaTemp.end());

	// -Update dt such that abs(Tp(n+1)-Tp(n)) <= Threshold
	double Threshold = 10;

	if (DeltaTemp_max > 0)
	{
		dt = dt_old * Threshold / DeltaTemp_max;
		dt = max(0.9 * dt_old, min(1.1 * dt_old, dt)); // Limit to only 10% change in dt
	}

	// -Adjust the remaining number of cycles according to the new dt
	dt = min(dt, globalTimeStep);  //in case of globalTimeStep>dt
	m_f = round(globalTimeStep / dt);
	if ((m_f * dt) != globalTimeStep) {    //in case of dt*m_f not equal globalTimeStep
		m_f = m_f + 1;
		dt = globalTimeStep / m_f;
	}
}


//Access

bool particle_1D::getState()
{
	return state;
}

double particle_1D::getNcells()
{
	return nx;
}

double particle_1D::getDelta()
{
	return delta;
}

double particle_1D::getMass()
{
	return mass_p;
}

double particle_1D::getVol()
{
	return volume;
}

double particle_1D::getMLR()
{
	return MLR;
}

double particle_1D::getGFRR()
{
	return GFRR;
}

double particle_1D::getCharHRR()
{
	return charHRR;
}

double particle_1D::getHconv()
{
	return h_conv;
}

std::vector<double> particle_1D::getCellCenters()
{
	return xCenter;
}

std::vector<double> particle_1D::getT()
{
	return Temp;
}

std::vector<double> particle_1D::getXm()
{
	return x_m;		
}

std::vector<double> particle_1D::getXvs()
{
	return x_vs;
}
