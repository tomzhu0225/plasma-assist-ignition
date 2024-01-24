function ctb = compute_third_body_concentration(Ck)

global n_reaction;
global n_species;
global third_body_reactions;
global third_body_coeff;
global third_body_species;
global species_name;

ctb = zeros(n_reaction,1);

for i=1:n_reaction
    if third_body_reactions(i) == 1
        for k=1:size(third_body_species,2)
            for j =1:n_species
                if (third_body_species(i,k) == species_name(j))
                   ctb(i) = ctb(i) + third_body_coeff(i,k) * Ck(j);
                end
            end
        end
    else
        ctb(i) = 1.0;
    end
end
end