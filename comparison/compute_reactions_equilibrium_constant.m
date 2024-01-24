function K_eq =  compute_reactions_equilibrium_constant(T)

global Patm;
global R;
global species_molar_mass;
global n_reaction;
global n_species;
global reactants_stochio;
global products_stochio;


K_eq = zeros(n_reaction,1);
delta_Z = zeros(n_reaction,1);
delta_coeff = zeros(n_reaction,1);

% compute species mass enthalpy
hkm = compute_species_mass_enthalpy(T).*species_molar_mass;

% compute species mass entropy
skm = compute_species_mass_entropy(T).*species_molar_mass;

for i=1:n_reaction
    for k=1:n_species
       delta_Z(i) = delta_Z(i) + (products_stochio(i,k)- reactants_stochio(i,k))*( (skm(k)/R) - (hkm(k)/(R*T)) );
       delta_coeff(i) = delta_coeff(i) + (products_stochio(i,k)- reactants_stochio(i,k));
    end
    K_eq(i) = (((Patm/(R*T))*10^(-6))^delta_coeff(i)) * exp(delta_Z(i)); 
end
end