function [sk] = compute_species_mass_entropy(T) 

global n_species;
global species_molar_mass;
global thermo_coeff;
global R ;
sk = zeros(n_species,1);

for i = 1:n_species
% compute molar Cpk
    for j = 1:4
        if (T < 1000)
           sk(i) = sk(i) + thermo_coeff(j+1,i) * (T^(j))/j;
        else
           sk(i) = sk(i) + thermo_coeff(j+8,i) * (T^(j))/j; 
        end
    end
    if (T < 1000)
        sk(i) = sk(i) + thermo_coeff(7,i) + thermo_coeff(1,i)* log(T) ;
    else
        sk(i) = sk(i) + thermo_coeff(14,i) + thermo_coeff(8,i)* log(T); 
    end
% compute mass sk
sk(i) = sk(i)/species_molar_mass(i); 
end
sk(:) = R*sk(:);
end