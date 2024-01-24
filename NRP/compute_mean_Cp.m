function cp = compute_mean_Cp(T,Yk)
global n_species;


Cpk = zeros(1:n_species);

Cpk = compute_species_mass_Cp(T);
cp =0.0;
for i=1:n_species
    cp = cp + Cpk(i)*Yk(i);
end
end