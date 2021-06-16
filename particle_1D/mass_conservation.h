//solve mass conservation
void particle_1D::mass_conservation(double x_O2_g)
{
	double xm_times_dV=0.0;
	double xvs_times_dV = 0.0;
	double xc_times_dV = 0.0;
	double store = 0.0;

	for (int i = 0; i < nx_old;i++)
	{
		xm_times_dV = (x_m_old[i] * dV_old[i]) * exp(-dt * R1.get_A() * exp(-R1.get_Ta() / Temp_old[i]));
		xvs_times_dV = (x_vs_old[i] * dV_old[i]) / (1 + dt * R2.get_A() * exp(-R2.get_Ta() / Temp_old[i]));

		if (eta_c != 0)
		{
			store = xvs_times_dV * dt * R2.get_A() * exp(-R2.get_Ta() / Temp_old[i])
				* (eta_c * rho_vs / rho_c);

			xc_times_dV = (x_c_old[i] * dV_old[i] + store)
				/ (1 + dt * x_O2_g * R3.get_A() * exp(-R3.get_Ta() / Temp_old[i]));
		}
		else
		{
			xc_times_dV = 0.0;
		}

		dV[i] = xm_times_dV + xvs_times_dV + xc_times_dV;

		x_m[i] = xm_times_dV / dV[i];
		x_m[i] = max(0.0, min(1.0, x_m[i])); 	// Safety: enforce 0 <= x_m <= 1

		x_vs[i] = xvs_times_dV / dV[i];
		x_vs[i] = max(0.0, min(1.0, x_vs[i]));	// Safety: enforce 0 <= x_vs <= 1

		x_c[i] = 1.0 - x_m[i] - x_vs[i];
		x_c[i] = max(0.0, min(1.0, x_c[i]));	// Safety: enforce 0 <= x_c <= 1

	}
}