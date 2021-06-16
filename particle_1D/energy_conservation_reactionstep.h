// Energy conservation
// - Sub-step associated with energy release due to pyrolysis/char oxidation
vector<double> particle_1D::energy_conservation_reactionstep(double x_O2_g)
{
	vector<double> Temp_return;
	Temp_return.assign(nx_old, 0.0);

	kp.assign(nx_old, 0.0);
	rho_times_cp.assign(nx_old, 0.0);
	Qdotp.assign(nx_old, 0.0);

	double dQdotpdT, RHS1_old, dRHS1dT_old, ratio;

	for (int i = 0; i < nx_old; i++)
	{
		kp[i] = k_m * x_m_old[i]
			  + k_vs * x_vs_old[i]
			  + k_c * x_c_old[i];

		rho_times_cp[i] = rho_m * c_m * x_m_old[i]
						+ rho_vs * c_vs * x_vs_old[i]
						+ rho_c * c_c * x_c_old[i];

		Qdotp[i] = rho_m * x_m_old[i] * R1.get_A() * exp(-R1.get_Ta() / Temp_old[i]) 
					* R1.get_deltaH()
				 + rho_vs * x_vs_old[i] * R2.get_A() * exp(-R2.get_Ta() / Temp_old[i]) 
					* R2.get_deltaH() * (1 - eta_c)
				 + rho_c * x_c_old[i] * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp_old[i]) 
					* R3.get_deltaH();

        dQdotpdT = rho_m * x_m_old[i] * R1.get_A() * exp(-R1.get_Ta() / Temp_old[i]) 
        			* R1.get_deltaH() * (R1.get_Ta()/pow(Temp_old[i],2.0))
        		 + rho_vs * x_vs_old[i] * R2.get_A() * exp(-R2.get_Ta() / Temp_old[i]) 
        		 	* R2.get_deltaH() * (1 - eta_c) * (R2.get_Ta()/pow(Temp_old[i],2.0))
        		 + rho_c * x_c_old[i] * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp_old[i]) 
        		 	* R3.get_deltaH() * (R3.get_Ta()/pow(Temp_old[i],2.0));

       RHS1_old    = Qdotp[i]/rho_times_cp[i];
       dRHS1dT_old = dQdotpdT/rho_times_cp[i];

        if( abs(dRHS1dT_old) > 0 )
        {
	       ratio   = RHS1_old/dRHS1dT_old;

	       Temp_return[i] = Temp0[i]
	                	  + (ratio+Temp0[i]-Temp_old[i])*(exp(0.5*dRHS1dT_old*dt)-1);
        }
        else // Case dRHS1dT_old = 0
        {
        	Temp_return[i] = Temp0[i] + 0.5*RHS1_old*dt;
        }
	}

	return Temp_return;
}