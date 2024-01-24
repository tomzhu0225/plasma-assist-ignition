function [h] = compute_mean_enthalpy(Yk,T) 

global n_species;
%global T0;

hk = compute_species_mass_enthalpy(T);

h = 0.0;

for i = 1:n_species
h = h + Yk(i) * hk(i); 
end
end