function M = compute_mean_molar_mass_from_mass_fractions(Yk)
global n_species;
global species_molar_mass;

M = 0.0;
for i = 1:n_species
    M = M + (Yk(i)/species_molar_mass(i));
end
M = 1/M;
end