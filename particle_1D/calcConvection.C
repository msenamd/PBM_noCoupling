#include "calcConvection.h"

double get_h(double T_g , double T_surf , double u_g , double delta , std::string geom) {

double Pr, nu_0, k_0, T_film , nu_g, k_g;
double D_eff, Re_D, Nu_D , h;
double nu_s, rho_g, rho_s, mu_g, mu_s;
const double pi=3.14159265358979323846;

Pr          = 0.7;      // Air Prandtl number
nu_0        = 1.6e-5;   // Air kinematic viscosity at ambient temp. [m2/s]
k_0         = 0.026;    // Air thermal conductivity at ambient temp.

T_film      = 0.5*(T_surf+T_g);         // Film temperature [K]

nu_g        = nu_0*pow((T_film/300.0),1.76); // Kin. visc. at film temp.
k_g         = k_0*pow((T_film/300.0),0.76);  // Conductivity at film temp.

D_eff       = 0.0;
Re_D        = 0.0;
Nu_D        = 0.0;


if(geom=="rectangle")   //rectangle
{
   D_eff       = (2*delta)*2/sqrt(pi); // Effective diameter of particle [m]
   Re_D        = u_g*D_eff/nu_g;       // Reynolds number at film temp.
   // Nusselt number at film temp. , Ref: Churchill and Bernstein
   Nu_D        = 0.3 + (0.62*pow(Re_D,0.5)*pow(Pr,0.3333))/pow((1+pow(0.4/Pr,2.0/3.0)),0.25)
                     * pow((1+pow(Re_D/282000,5.0/8.0)),4.0/5.0) ;
}
else if(geom=="cylinder")   //cylinder
{
   D_eff       = 2*delta;              // Effective diameter of particle [m]
   Re_D        = u_g*D_eff/nu_g;       // Reynolds number at film temp.
   // Nusselt number at film temp. , Ref: Churchill and Bernstein
   Nu_D        = 0.3 + (0.62*pow(Re_D,0.5)*pow(Pr,0.3333))/pow((1+pow(0.4/Pr,2.0/3.0)),0.25)
                     * pow((1+pow(Re_D/282000,5.0/8.0)),4.0/5.0) ;
}
else if(geom=="sphere")  //sphere
{
   D_eff       = 2*delta;                             // Effective diameter of particle [m]
   k_g         = k_0*pow((T_g/300.0),(0.76));         // Conductivity at gas temp.
   nu_g        = nu_0*pow((T_g/300.0),(1.76));        // Kin. visc. at gas temp.
   nu_s        = nu_0*pow(T_surf/300.0 , 1.76);  // Kin. visc. at surface temp.
   rho_g       = 101325.0/(287.0*T_g);
   rho_s       = 101325.0/(287.0*T_surf);
   mu_g        = rho_g*nu_g;
   mu_s        = rho_s*nu_s;
   Re_D        = u_g*D_eff/nu_g;                      // Calculated at gas temperature
   Nu_D        = 2+(0.4*pow(Re_D,0.5)+0.06*pow(Re_D, 2.0/3.0))*pow(Pr,0.4)*pow(mu_g/mu_s,0.25);
}

h    = Nu_D*k_g/D_eff;   // Convective heat transfer coef. [W/m2/K]

return h;

}
