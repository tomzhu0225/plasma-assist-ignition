% Main script

% Declare global Variables
global R ;
global n_species;
global species_name;
global P;
global Patm;
global Yk_FG;
global time;
global tau
global T_in
global Yk_in

% Set global varibales values     !!!!!! DON'T CHANGE THIS BLOC
R = 8.314;
Patm = 101325;
tau=10e-3
T_in=300

% read species thermal properties        !!!!!! DON'T CHANGE THIS BLOC
file_name_thermo = 'H2_thermo.txt';
read_thermo(file_name_thermo);
file_name_chemistry = 'H2_chemistry.txt';
read_chemistry(file_name_chemistry);

% Declare variables                      !!!!!! DON'T CHANGE THIS BLOC
Yk_FG = zeros(n_species,1);
z0 = zeros(n_species+1,1);
result = zeros(1,n_species+2);

% Set time control
final_time = 50e-03;                % Simulation total time  
time = 0.0;                        % Simulation initial time
dt = 10e-05;                       % time step

equiv=0.35
gamma=3.76
% examples on how to use the functions can be found in " examples.m" script
Yk_FG(4)= 0.002*equiv/(0.002*equiv+0.5*(0.032+gamma*0.028));  
Yk_FG(2)= 0.032*0.5/(0.002*equiv+0.5*(0.032+gamma*0.028))

Yk_FG(6)= 0;
Yk_FG(9)=1-Yk_FG(2)-Yk_FG(4);


Yk_in=Yk_FG

% Run the Simulation 



% Set case variables
P = 101325 ;    % Pressure
T_FG =  1000  % initial temperature

% intialization of z0 (T, Yk) 
z0(1) = T_FG
z0(2:end) = Yk_FG
% The time, temperature and species are stored in the matrix result
% you can use this matrix to plot the different variables
result = [time transpose(z0)];

while (time < final_time)
    
tspan = [0 dt];
% Call the LSODI solver

options = odeset('RelTol',1e-6,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15 ...
    1e-15 1e-15 1e-15 1e-15]);

% z0 is the intial solution at t
% z(size(z,1),:) is the solution at t+dt

[t,z] = ode15s(@ode,tspan,z0,options);

result = [result; t+time z];    % Store the result ( Time, T , Yk )

z0 = z(size(z,1),:);            % Update   z0
time = time + dt;               % Update time

end


time = result(:, 1);        % First column is time
temperature = result(:, 2); % Second column is temperature

% Main species columns
H2 = result(:, 6);
O2 = result(:, 4);
H2O = result(:, 8);
N2 = result(:, 11);

% Minor species columns
O = result(:, 3);
H = result(:, 5);
OH = result(:, 7);
HO2 = result(:, 9 );
H2O2 = result(:,10);                     

% Plot Time vs Temperature
figure;
plot(time, temperature, 'LineWidth', 2);
title('Time vs Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');
grid on;

% Plot Time vs Main Species
figure;
plot(time, O2, 'b-', 'LineWidth', 2); hold on;
plot(time, H2, 'r-', 'LineWidth', 2); 
plot(time, H2O, 'g-', 'LineWidth', 2);
plot(time, N2, 'k-', 'LineWidth', 2);
legend({'O2', 'H2', 'H2O', 'N2'});
title('Time vs Main Species');
xlabel('Time (s)');
ylabel('Mass Fraction');
grid on;

% Plot Time vs Minor Species
figure;
plot(time, O, 'b-', 'LineWidth', 2); hold on;
plot(time, H, 'r-', 'LineWidth', 2);
plot(time, OH, 'g-', 'LineWidth', 2);
plot(time, HO2, 'k-', 'LineWidth', 2);
plot(time, H2O2, 'm-', 'LineWidth', 2);
legend({'O', 'H', 'OH', 'HO2', 'H2O2'});
title('Time vs Minor Species');
xlabel('Time (s)');
ylabel('Mass Fraction');
grid on;



