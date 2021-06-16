// Energy conservation
// - Sub-step associated with heat conduction plus gas-to-solid heat
//   transfer at the exposed surface
// - Implicit time integration (Crank-Nicolson)
//   (implicit treatment of both diffusion and of gas-to-solid heat transfer)
// - Formulation as a tri-diagonal matrix system
void particle_1D::energy_conservation_diffusionstep(double T_g, double G, double x_O2_g)
{
	// constructing vectors
	S_pos.assign(nx_old, 0.0);
	S_neg.assign(nx_old, 0.0);

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

	//Constructing Temperature matrix (filling vectors A,B,C,D)
	int i;

    // innermost surface boundary cell
    i = 0;
    FO_L = 0;
    set_FO_R(i, dt);

	vectA[i] = -0.5*FO_L;
	vectB[i] = 1 + 0.5*FO_L + 0.5*FO_R;
	vectC[i] = -0.5*FO_R;
	vectD[i] = (1-0.5*FO_R)*Temp_old[i] + 0.5*FO_R*Temp_old[i+1];

    // exposed surface boundary cell
    i = nx_old - 1;
    FO_R = 0;
    set_FO_L(i, dt);

    vectA[i] = -0.5*FO_L;
    vectB[i] = 1 + 0.5*FO_L + 0.5*FO_R;
    vectC[i] = -0.5*FO_R;
    vectD[i] = (1 - 0.5*FO_L) * Temp_old[i] + 0.5*FO_L*Temp_old[i-1];

    kp_surf  = kp.back();

    eps_surf = e_m * x_m.back()
             + e_vs * x_vs.back()
             + e_c * x_c.back();

	dx_surf = xRight.back() - xCenter.back();

	double h_rad = eps_surf * sigma * pow(Tsurf_old, 3.0);

	double Bi = (h_conv + h_rad) * dx_surf / kp_surf;

	vectB[i] = vectB[i] + ((h_conv + h_rad) / (1+Bi))
            	* 0.5*dt/rho_times_cp[i] * S_pos[i]/dV[i];

	vectD[i] = vectD[i] + ( q_surf_old + ((eps_surf*G + h_conv*T_g)/(1+Bi)) )
            	* 0.5*dt/rho_times_cp[i] * S_pos[i]/dV[i];

    // interior cells
    for (i = 1; i < (nx_old - 1); i++)
    {
        set_FO_L(i, dt);
        set_FO_R(i, dt);

        vectA[i] = -0.5*FO_L;
       	vectB[i] = 1 + 0.5*FO_L + 0.5*FO_R;
        vectC[i] = -0.5*FO_R;
        vectD[i] = (1 - 0.5*FO_L - 0.5*FO_R)*Temp_old[i]
          			+ 0.5*FO_L*Temp_old[i-1] + 0.5*FO_R*Temp_old[i+1];
    }
}