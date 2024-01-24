function read_thermo(file_name)

global n_elements;
global elements_name;
global elements_molar_mass;
global n_species;
global species_name;
global species_molar_mass;
global thermo_coeff;
global species_composition;


fid=fopen(file_name,'r') ;
fscanf(fid, '%s',1);

% Read elements 

n_elements = fscanf(fid, '%d',1) ;
elements_name = strings(1,n_elements);
elements_molar_mass = zeros(n_elements);
dummy_str = fscanf(fid, '%s',n_elements);
for k = 1:n_elements
  elements_name(k) = dummy_str(k);
end
elements_molar_mass = fscanf(fid, '%f',n_elements)*1E-3;

% Read species

fscanf(fid, '%s',1);
n_species = fscanf(fid, '%d',1) ;
species_name = strings(1,n_species);
species_table = zeros(n_species,6);
species_composition = table('Size',[n_species 3],'VariableTypes',{'double','double','double'},'VariableNames',{'H','O','N'});
thermo_coeff = zeros(14,n_species);
species_molar_mass = zeros(n_species,1);

for i = 1:n_species
    species_name(i) = fscanf(fid, '%s',1);
    species_table(i,:) = fscanf(fid, '%s%f%s%f%s%f',6);
    species_molar_mass (i) = fscanf(fid, '%f',1) * 1E-3;
    thermo_coeff(:,i) = fscanf(fid, '%f',14);
end
species_composition.Properties.RowNames = species_name;
species_composition.H = species_table(:,2);
species_composition.O = species_table(:,4);
species_composition.N = species_table(:,6);
fclose(fid);

% Clear variables
clearvars species_table;
clearvars i;
clearvars k;
clearvars fid;
clearvars dummy_str

end

