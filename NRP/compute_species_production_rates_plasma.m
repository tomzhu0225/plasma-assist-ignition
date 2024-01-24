function wdot =  compute_species_production_rates_plasma(Ep,Yk,Yk_in,hk)

global species_molar_mass;
eta=0.63;

wdot=zeros(length(species_molar_mass),1);
wdot(1)=eta*Yk(2)/Yk_in(2)*(Ep/(species_molar_mass(1)*(hk(1)-hk(2))));
wdot(2)=-species_molar_mass(1)/species_molar_mass(2) * wdot(1);

end