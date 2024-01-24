% Main script

% Declare global Variables
global R ;
global n_species;
global species_name;
global P;
global Patm;
global Yk_FG;
global time;


% Set global varibales values     !!!!!! DON'T CHANGE THIS BLOC
R = 8.314;
Patm = 101325;

% read species thermal properties        !!!!!! DON'T CHANGE THIS BLOC
file_name_thermo = 'H2_thermo.txt';
read_thermo(file_name_thermo);
file_name_chemistry = 'H2_chemistry.txt';
read_chemistry(file_name_chemistry);

% Declare variables                      !!!!!! DON'T CHANGE THIS BLOC


equivdelay_temp=[]
for T_FG =  1000:150:1300
T_FG
equivdelay=[]
for equiv=0.5:0.1:2
equiv
Yk_FG = zeros(n_species,1);
z0 = zeros(n_species+1,1);

% Set time control
final_time = 4e-03;                % Simulation total time  
time = 0.0;                        % Simulation initial time
dt = 2e-06;                       % time step
gamma=3.76;


% examples on how to use the functions can be found in " examples.m" script
Yk_FG(4)= 0.002*equiv/(0.002*equiv+0.5*(0.032+gamma*0.028));  
Yk_FG(2)= 0.032*0.5/(0.002*equiv+0.5*(0.032+gamma*0.028));

Yk_FG(6)= 0;
Yk_FG(9)=1-Yk_FG(2)-Yk_FG(4);
% Run the Simulation 

% TO COMPLETE

% Set case variables
P = 101325 ;    % Pressure
 % initial temperature

% intialization of z0 (T, Yk) 
z0(1) = T_FG;
z0(2:end) = Yk_FG;
% The time, temperature and species are stored in the matrix result
% you can use this matrix to plot the different variables
result = [time transpose(z0)];
wH2=1;

while 1

tspan = [0 dt];
% Call the LSODI solver

options = odeset('RelTol',1e-6,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15 ...
    1e-15 1e-15 1e-15 1e-15]);

% z0 is the intial solution at t
% z(size(z,1),:) is the solution at t+dt

[t,z] = ode15s(@ode,tspan,z0,options);

% Find elements smaller than 'a'
vector= z(:,5);
smallerThanA = vector(vector<0.99*Yk_FG(4));
if isempty(smallerThanA);
    fprintf('No value smaller than %d found in the vector.\n');
else
    % Find the value closest to 'a' and its index in the original vector
    [~, idxInSmaller] = min(abs(smallerThanA - Yk_FG(4)));
    closestValue = smallerThanA(idxInSmaller);
    originalIndex = find(z(:,5) == closestValue, 1);
    t_delay=t+time;
    t_delay=t_delay(originalIndex);
    break
end

result = [result; t+time z];    % Store the result ( Time, T , Yk )

z0 = z(size(z,1),:);            % Update   z0
time = time + dt;               % Update time

end

equivdelay=[equivdelay;equiv t_delay];
end
equivdelay_temp=[equivdelay_temp equivdelay];
end
x = equivdelay_temp(:,1);  % Equivalence ratio
y1 = equivdelay_temp(:,2); % Autoignition delay at 1000K
y2 = equivdelay_temp(:,4); % Autoignition delay at 1150K
y3 = equivdelay_temp(:,6); % Autoignition delay at 1300K

% Create the plot
plot(x, y1, 'b-', x, y2, 'r-', x, y3, 'g-');
xlabel('Equivalence Ratio');
ylabel('Autoignition Delay (s)');
title('Autoignition Delay vs. Equivalence Ratio at Different Temperatures');
legend('1000K', '1150K', '1300K');

% Adjust additional plot settings as necessary
grid on;


