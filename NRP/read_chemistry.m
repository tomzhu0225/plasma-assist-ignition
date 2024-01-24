function read_chemistry(file_name)
global n_reaction;
global A_cst;
global T_power;
global AE;
global n_species;
global species_name;
global reactants_stochio;
global products_stochio;
global third_body_species;
global third_body_coeff;
global third_body_reactions;

%n_species = 9;

%file_name = 'H2_chemistry.txt';
%file_name = 'warnatz.txt';
fid=fopen(file_name,'r') ;
line_1 = convertCharsToStrings(fgetl(fid));
line_1 = split(line_1);
n_reaction = str2num(line_1(2));
A_cst = zeros(n_reaction,1);
T_power = zeros(n_reaction,1);
AE = zeros(n_reaction,1);
reactants_stochio = zeros(n_reaction,n_species);
products_stochio = zeros(n_reaction,n_species);
third_body_reactions = zeros(n_reaction,1);
third_body_species = strings(n_reaction,4);



for i = 1:n_reaction
    line = convertCharsToStrings(fgetl(fid));
    line = split(line);
    
    % reaction parameters
    A_cst(i) = str2double(line(2));
    T_power(i) = str2double(line(3));
    AE(i) = str2double(line(4));
    
    % reactants and products
    A = split(line(1),'<=>');
    B = split(A(1),'+');
    reactants(i,1:size(B,1)) = B;
    B = split(A(2),'+');
    products(i,1:size(B,1)) = B;
    
    % species  stochio coeff
    % Reactants
    for k = 1:size(reactants,2)
        react = reactants(i,k);
        react = convertStringsToChars(react);
        if ~isempty(react)
            if (react(1) == 'H' || react(1) == 'O' || react(1) == 'N' || react(1) == 'M')
                for j = 1:n_species
                    if (react == species_name(j))
                        reactants_stochio(i,j) = 1.0;
                    end
                end
                if (react == 'M')
                    third_body_reactions(i) = 1;
                    line_M = convertCharsToStrings(fgetl(fid));
                    line_M = split(line_M);
                    C = split(line_M,'/');
                    third_body_species(i,1:4)= C(1:size(C,1),1);
                    third_body_coeff(i,1:size(C,1))= C(1:size(C,1),2);
                end
            else
                for j = 1:n_species
                    if (react(2:length(react)) == species_name(j))
                        reactants_stochio(i,j) = str2num(react(1)) ;
                    end
                end
            end
        end
    end
    
    % Products
    for k = 1:size(products,2)
        prod = products(i,k);
        prod = convertStringsToChars(prod);
        if ~isempty(prod)
            if (prod(1) == 'H' || prod(1) == 'O' || prod(1) == 'N' || prod(1) == 'M')
                for j = 1:n_species
                    if (prod == species_name(j))
                        products_stochio(i,j) = 1.0;
                    end
                end
            else
                for j = 1:n_species
                    if (prod(2:length(prod)) == species_name(j))
                        products_stochio(i,j) = 1* str2num(prod(1)) ;
                    end
                end
            end
        end
    end
    
end
clearvars A;
clearvars B;
clearvars reactants;
clearvars products;
clearvars react;
clearvars prod;
clearvars i;
clearvars j;
clearvars k;
clearvars fid;
clearvars line_1;
clearvars line;
clearvars line_M;
clearvars C;
end