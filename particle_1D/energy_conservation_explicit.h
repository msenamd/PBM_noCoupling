// Energy conservation
// - Explicit time integration
vector<double> particle_1D::energy_conservation_explicit(double x_O2_g)
{
	vector<double>Temp_explicit;
	Temp_explicit.assign(nx_old, 0.0);

	S_pos.assign(nx_old, 0.0);
	S_neg.assign(nx_old, 0.0);

	kp.assign(nx_old, 0.0);
	rho_times_cp.assign(nx_old, 0.0);
	Qdotp.assign(nx_old, 0.0);

    if (shape == "rectangle")
	{
        S_pos[0] = areaRectangle;
        S_neg[0] = areaRectangle;
		for (int i = 1; i < nx_old; i++)
		{
            S_pos[i] = areaRectangle;
            S_neg[i] = areaRectangle;
		}
	}
    else if (shape == "cylinder")
	{
        S_pos[0] = 2.0 * pi * xRight[0] * lengthCylinder;
		S_neg[0] = 0.0;
		for (int i = 1; i < nx_old; i++)
		{
            S_pos[i] = 2.0 * pi * xRight[i] * lengthCylinder;
            S_neg[i] = 2.0 * pi * xRight[i - 1] * lengthCylinder;
		}
	}
    else if (shape == "sphere")
	{
		S_pos[0] = 4.0 * pi * pow(xRight[0], 2.0);
		S_neg[0] = 0.0;
		for (int i = 1; i < nx_old; i++)
		{
			S_pos[i] = 4.0 * pi * pow(xRight[i], 2.0);
			S_neg[i] = 4.0 * pi * pow(xRight[i - 1], 2.0);
		}
	}

	for (int i = 0; i < nx_old; i++)
	{
        kp[i] = k_m * x_m_old[i]
              + k_vs * x_vs_old[i]
              + k_c * x_c_old[i];

        rho_times_cp[i] = rho_m * c_m * x_m_old[i]
            			+ rho_vs * c_vs * x_vs_old[i]
            			+ rho_c * c_c * x_c_old[i];

        Qdotp[i] = rho_m * x_m_old[i] * R1.get_A() * exp(-R1.get_Ta() / Temp_old[i]) * R1.get_deltaH()
            	 + rho_vs * x_vs_old[i] * R2.get_A() * exp(-R2.get_Ta() / Temp_old[i]) * R2.get_deltaH() * (1 - eta_c)
            	 + rho_c * x_c_old[i] * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp_old[i]) * R3.get_deltaH();

	}

	int i;

	// innermost surface boundary cell
    i = 0;
    set_FO_R(i, dt_old);
    Temp_explicit[i] = (1 - FO_R) * Temp_old[i] + FO_R * Temp_old[i + 1]
            		 + Qdotp[i] * dt_old / rho_times_cp[i];

    // exposed surface boundary cell
    i = nx_old - 1;
    set_FO_L(i, dt_old);
    Temp_explicit[i] = FO_L * Temp_old[i - 1] + (1 - FO_L) * Temp_old[i]
            		 + dt_old / rho_times_cp[i] * Qdotp[i]
            		 + dt_old / rho_times_cp[i] * S_pos[i] / dV[i] * q_surf_old;

    // interior cells
    for (i = 1; i < (nx_old - 1); i++)
    {
            set_FO_L(i, dt_old);
            set_FO_R(i, dt_old);

            Temp_explicit[i] = FO_L * Temp_old[i - 1] + (1 - FO_L - FO_R) * Temp_old[i]
                			 + FO_R * Temp_old[i + 1] + dt_old / rho_times_cp[i] * Qdotp[i];
    }

    return Temp_explicit;
}
