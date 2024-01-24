function hk = compute_species_mass_enthalpy(T) 

global n_species;
global species_molar_mass;
global thermo_coeff;
global R ;
hk = zeros(n_species,1);

for i = 1:n_species
% compute molar Cpk
    for j = 1:5
        if (T < 1000)
           hk(i) = hk(i) + thermo_coeff(j,i) * (T^(j))/j;
        else
           hk(i) = hk(i) + thermo_coeff(j+7,i) * (T^(j))/j; 
        end
    end
    if (T < 1000)
        hk(i) = hk(i) + thermo_coeff(6,i);
    else
        hk(i) = hk(i) + thermo_coeff(13,i); 
    end
% compute mass hk
hk(i) = hk(i)/species_molar_mass(i);
end
hk(:) = R*hk(:);
end