function dz = ode(t,z,flag)
global n_species;


% set variables
temperature = z(1);
Yk(1:n_species) = z(2:n_species+1);

dz = auto_ignition_cst_pressure(temperature,Yk,flag);



end 