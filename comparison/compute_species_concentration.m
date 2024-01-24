function Ck = compute_species_concentration(rho,Yk)
global n_species;
global species_molar_mass;

Ck = zeros(1,n_species);
for i=1:n_species
    Ck(i) = rho * Yk(i)/species_molar_mass(i);
end
end