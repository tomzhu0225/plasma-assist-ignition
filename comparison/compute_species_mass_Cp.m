function [Cpk] = compute_species_mass_Cp(T) 

global n_species;
global thermo_coeff;
global species_molar_mass;
global R ;
Cpk = zeros(n_species,1);

for i = 1:n_species
% compute molar Cpk
    for j = 1:5
        if (T < 1000)
           Cpk(i) = Cpk(i) + thermo_coeff(j,i) * T^(j-1);
        else
           Cpk(i) = Cpk(i) + thermo_coeff(j+7,i) * T^(j-1); 
        end
    end
% compute mass Cpk
Cpk(i) = Cpk(i)/ species_molar_mass(i);
end
Cpk(:) = R*Cpk(:);

end