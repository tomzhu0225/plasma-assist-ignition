function dz = auto_ignition_cst_pressure(temperature,Yk,pulse_flag)

global n_species;
global P;
global Ck;
global species_molar_mass;
global tau;
global T_in;
global Yk_in;
global Ed
global Vd
global td


Ep=Ed/(Vd*td)

dz= zeros(n_species+1,1);
rho = compute_density_from_pressure(P,temperature,Yk);

Ck=rho*Yk'./species_molar_mass;

wdot =  compute_species_production_rates(temperature,Ck);
hk = compute_species_mass_enthalpy(temperature);
hk_in = compute_species_mass_enthalpy(T_in);
wdotp=compute_species_production_rates_plasma(Ep,Yk,Yk_in,hk);
cp = compute_mean_Cp(temperature,Yk);
dz(2:10)= species_molar_mass .* wdot/rho+ (Yk_in-Yk')/tau+ pulse_flag*species_molar_mass .* wdotp/rho;
dz(1)=-sum(hk.*species_molar_mass.*wdot)/(rho*cp)+(1/(tau*cp))*(sum(hk_in.*Yk_in)-sum(Yk_in.*hk))+pulse_flag*(Ep/(cp*rho)-sum(hk.*species_molar_mass.*wdotp)/(rho*cp));


% To complete ...

end