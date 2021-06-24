// Remeshing steps
// Design objective:
// - Maintain uniform grid with spacing close to initial value

void particle_1D::remeshing()
{

// Remeshing step 1: update nx, dx
	xCenter_old = xCenter;
	xRight_old 	= xRight;
	dx_old 		= dx;
    dV_old 		= dV;

	// Update the number of cells
	nx = round(xRight_old.back() / dx_i); 	// Update the number of cells
	nx = min(nx, nx_old);   				// nx <= nx_old
	nx = max(nx, 5);						// Keep minimum of 5 grid cells

	dx.assign(nx, ((xRight_old.back() / nx)));

// Remeshing, step 2: update computational grid
	
	xRight 	= get<0>(setMesh(nx, dx, shape, areaRectangle, lengthCylinder));
	xCenter = get<1>(setMesh(nx, dx, shape, areaRectangle, lengthCylinder));
	dV 		= get<2>(setMesh(nx, dx, shape, areaRectangle, lengthCylinder));

	volume 	= accumulate(dV.begin(), dV.end(), 0.0);

// Remeshing, step 3: interpolate solution on new computational grid
	Temp_old 	= Temp;
	x_m_old 	= x_m;
	x_vs_old 	= x_vs;
	x_c_old 	= x_c;

	interpolateOnNewMesh(nx);	

}


// Look for q that conserves int_{x=0}^{x=volume}(q dV)
// q = x_m, x_vs, x_c, (rho_cp_T)
void particle_1D::interpolateOnNewMesh(int nx_new)
{

	vector<double> rho_cp_T_old, rho_cp_T_new;
	vector<double> x_m0, x_vs0, x_c0, rho_cp_T0;

	Temp.assign(nx_new,0.0);
	x_m.assign(nx_new, 0.0);
	x_vs.assign(nx_new, 0.0);
	x_c.assign(nx_new, 0.0);

	x_m0.assign(nx_new, 0.0);
	x_c0.assign(nx_new, 0.0);
	x_vs0.assign(nx_new, 0.0);
	rho_cp_T0.assign(nx_new, 0.0);

	rho_cp_T_old.assign(nx_old, 0.0);
	rho_cp_T_new.assign(nx_new, 0.0);

	for (int i = 0; i < nx_old; i++)
	{
	    rho_cp_T_old[i] = (rho_m*c_m*x_m_old[i] + rho_vs*c_vs*x_vs_old[i] 
	    								+ rho_c*c_c*x_c_old[i]) * Temp_old[i];		
	}

	// First estimate of q based on linear interpolation: q0
	// q0 = x_m0, x_vs0, x_c0, (rho_cp_T0)
	int ii=-1;
	for (int i = 0; i < nx; i++)
	{
		ii = -1;
		for (int j = 0; j < (nx_old - 1); j++)
		{
			if (((xCenter_old[j] - xCenter[i]) <= 0.0) && ((xCenter[i] - xCenter_old[j+1]) < 0.0))
			{
				ii = j;
			}

		}
		if ((ii == -1) && (xCenter[i] < xCenter_old[0]))
		{
			ii = 0;
		}
		if ((ii == -1) && (xCenter_old[nx_old-1] <= xCenter[i]))
		{
			ii = nx_old - 2;
		}
		if (ii > -1)
		{
			rho_cp_T0[i] = rho_cp_T_old[ii]
				+ (rho_cp_T_old[ii + 1] - rho_cp_T_old[ii]) * (xCenter[i] - xCenter_old[ii])
				/ (xCenter_old[ii + 1] - xCenter_old[ii]);

			x_m0[i] = x_m_old[ii]
				+ (x_m_old[ii + 1] - x_m_old[ii]) * (xCenter[i] - xCenter_old[ii])
				/ (xCenter_old[ii + 1] - xCenter_old[ii]);

			x_vs0[i] = x_vs_old[ii]
				+ (x_vs_old[ii + 1] - x_vs_old[ii]) * (xCenter[i] - xCenter_old[ii])
				/ (xCenter_old[ii + 1] - xCenter_old[ii]);

			x_c0[i] = x_c_old[ii]
				+ (x_c_old[ii + 1] - x_c_old[ii]) * (xCenter[i] - xCenter_old[ii])
				/ (xCenter_old[ii + 1] - xCenter_old[ii]);
		}
		if (ii == -1)
		{
			cout << "*** Error in linear interpolation intial estimate ***" << endl;
			break;
		}
	}

	x_m = solveCostFunction(x_m0, x_m_old, nx_new, nx_old);
	checkBounds(x_m);

	x_vs = solveCostFunction(x_vs0, x_vs_old, nx_new, nx_old);
	checkBounds(x_vs);

	x_c = solveCostFunction(x_c0, x_c_old, nx_new, nx_old);
	checkBounds(x_c);

	rho_cp_T_new = solveCostFunction(rho_cp_T0, rho_cp_T_old, nx_new, nx_old);

	for (int i = 0; i < nx_new; i++)
	{
		Temp[i] = rho_cp_T_new[i]/ 
				  (rho_m*c_m*x_m[i] + rho_vs*c_vs*x_vs[i] + rho_c*c_c*x_c[i]);
	}
}


/* Find q that minimizes cost function:
   CF = W_I*(int_q - int_q_old)^2 ...
      + W_P*sum_{i=1}^{i=nx_new}(q(i)-qref(i))^2
   where int_q = sum_{i=1}^{i=nx_new}(q(i)*coef(i)*dV(i))
   and where W_I and W_P are weight coefficients
   Use Method of Conjugate Gradients  
*/ 
vector<double> particle_1D::solveCostFunction(vector<double> q0, vector<double> q_old, int nx_new, int nx_old)
{

	int iter_max = 10;
	double eps = 1e-6;

	// A dummy empty vector of size nx_new
	vector<double> zeroVector;
	zeroVector.assign(nx_new, 0.0);
	
	// Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
	double int_q_old = 0.0;
	for(int i=0; i<nx_old; i++)
	{
	    int_q_old += q_old[i]*dV_old[i];
	}

	// Scaling of CG coefficents
	double W_I      = 1e+6;
	double W_P      = 1e-3;
	double Q        = int_q_old/volume;
	Q        = max(Q,1e-10);
	W_I      = W_I*(1*0.01*0.0005/Q/volume/dV[1]);
	W_P      = W_P*(1/Q);

	// Reference solution
	vector<double> qref = q0;

	// Initial guess (taken as equal to reference solution)
	vector<double> q = zeroVector;
	vector<double> coef = zeroVector;

	for (int i=0; i<nx_new; i++)
	{
		if(qref[i] < 1e-10)
	    {
	        coef[i] = 0.0;
	    }
	    else
	    {
	    	coef[i] = 1.0;
	    }

		q[i] = q0[i] + coef[i]*1e-5*Q*(2*0.678-1); 	// Add noise
    	//q[i] = max(0.0,min(1.0,q[i]));            // Safety: 0 <= x_m <= 1
	}

	// Calculate initial residual
	int iter = 0;

	vector<double> res 	= zeroVector;
	vector<double> d 	= zeroVector;

	double int_q0 			= 0.0;
	double magnitude_res 	= 0.0;

	for (int i = 0; i < nx_new; i++)
	{
		int_q0 += q[i] * coef[i] * dV[i];
	}
	for (int i=0; i<nx_new; i++)
	{
		res[i] = -2*W_I*(int_q0 - int_q_old)*coef[i]*dV[i]
           		 -2*W_P*(q[i]-qref[i]);

    	d[i]   = res[i];

		magnitude_res += res[i]*res[i];
	}

	double magnitude_res0 = magnitude_res;

	double int_q = 0.0;

	if(magnitude_res0 <= 1e-14)
	{
		// Initial guess is the solution
	    int_q = int_q0;
	}
	else
	{
		vector<double> Ad = zeroVector;

		while((iter <= iter_max) && (magnitude_res/magnitude_res0 > pow(eps,2.0)) )
		{
			iter++;

			for (int i = 0; i < nx_new; i++)
			{
				Ad[i] = (2*W_P)*d[i];
		        for (int j = 0; j < nx_new; j++)
		        {
		            Ad[i] = Ad[i] + (2*W_I*coef[i]*dV[i]*coef[j]*dV[j])*d[j];
		        }	
			}

			double alpha = 0.0;
			for (int i = 0; i < nx_new; i++)
			{
				alpha = alpha + d[i] * Ad[i];
			}
			alpha = magnitude_res/alpha;

			int_q = 0;
			for (int i = 0; i < nx_new; i++)
			{
				q[i] = q[i] + alpha*d[i];

				int_q = int_q + q[i]*coef[i]*dV[i];
			}

        	for (int i = 0; i < nx_new; i++)
        	{
            	res[i] = -2*W_I*(int_q - int_q_old)*coef[i]*dV[i]
                    	 -2*W_P*(q[i]-qref[i]);
        	}

        	double magnitude_res_old = magnitude_res;
        	magnitude_res = 0;

        	for (int i = 0; i < nx_new; i++)
        	{
        		magnitude_res = magnitude_res + res[i]*res[i];
        	}

        	double beta = magnitude_res/magnitude_res_old;

        	for (int i = 0; i < nx_new; i++)
        	{
        		d[i] = res[i] + beta*d[i];
        	}
		}
		if (iter >= iter_max)
		{
			cout << "*** Error in remeshing cost function ***" << endl;
		}
	}

	return q;
}


//sets vector bounds to 0-1
void particle_1D::checkBounds(vector<double> q)
{
	for(unsigned int i=0; i< q.size(); i++)
	{
		q[i] = max(0.0, min(1.0, q[i]));
	}

}

