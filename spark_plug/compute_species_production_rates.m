function wdot =  compute_species_production_rates(T,Ck)

% Function to compute the species production rate wdot [mol/m3/s]

global n_reaction;
global A_cst;
global T_power;
global AE;
global R;
global n_species;
global reactants_stochio;
global products_stochio;
global species_molar_mass;
global ccf 
global ccr
global q
global ctb
global kf 
global kr
global keq
global wdot

Ck = Ck * 10^(-6);

% compute reactions rate constants
kf = zeros(n_reaction,1);
kr = zeros(n_reaction,1);

% compute equilibrium constants
keq =  compute_reactions_equilibrium_constant(T);

% compute foraward and reverse rate constants
% Arrhenius Law
for i=1:n_reaction
    kf(i) = A_cst(i) * (T^(T_power(i))) * exp(-AE(i)/(R*0.239*T));
    kr(i) = kf(i)/keq(i);
end
% compute reactions rate of progress
q = zeros(n_reaction,1);
ccf = zeros(n_reaction,1);
ccr = zeros(n_reaction,1);
for i = 1:n_reaction
    ccf(i) = 1.0;
    ccr(i) = 1.0;
end

% compute third bodies concentrations
ctb = compute_third_body_concentration(Ck);

% compute rate of progress
for i=1:n_reaction
    for k = 1:n_species
        ccf(i) = ccf(i) * Ck(k)^(reactants_stochio(i,k));
        ccr(i) = ccr(i) * Ck(k)^(products_stochio(i,k));
    end
    q(i) = ctb(i) *( (kf(i)*ccf(i)) - (kr(i)*ccr(i)));
end

% compute species production rate 
wdot = zeros(n_species,1);

for k = 1:n_species
    for i =1:n_reaction
       wdot(k) = wdot(k) + (products_stochio(i,k)- reactants_stochio(i,k)) * q(i);
    end
    wdot(k) = wdot(k) * 10^(6);
end

end