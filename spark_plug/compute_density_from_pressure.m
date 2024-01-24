function rho = compute_density_from_pressure(P,temperature,Yk)
global R;

M = compute_mean_molar_mass_from_mass_fractions(Yk);
rho = (P * M)/(temperature * R);
end 