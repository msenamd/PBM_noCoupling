void particle_1D::updateSurface(double T_g, double G)
{
	int iter = 0;
	int maxIter = 100;
	double error = 1.0;
	double T_old,T_new, h_rad, Bi;
		
	T_old = Tsurf_old;

	while ( (error > 1e-3) && (iter < maxIter) )
	{
		iter = iter + 1;
		h_rad = eps_surf * sigma * pow(T_old , 3.0);
		Bi = (h_conv + h_rad) * dx_surf / kp_surf;
		T_new = (Temp.back() + h_conv * T_g * dx_surf / kp_surf + (eps_surf * G * dx_surf / kp_surf)) / (1.0 + Bi);

		error = abs(T_new - T_old);
		T_old = T_new;
	}
	if (iter >= maxIter)
	{
		cout << "*** Error in calculating particle surface Temperature" <<
						", iterations exceed limit ***" << endl;		
		return;
	}

	Tsurf = T_new;
	
	q_surf = eps_surf * G - eps_surf * sigma * pow(Tsurf, 4.0)
					+ h_conv * (T_g - Tsurf);		
}



