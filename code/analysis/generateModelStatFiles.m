function generateModelStatFiles(model, out_dir)
%% generateModelStatsFile(model, out_dir)
% produces:
% > a data file 'model_stats.txt' with numbers for model features:
% -- reactions
% -- metabolites
% -- genes
% -- reactions with GPR
% -- metabolic reactions without GPR
% -- transport reactions
% -- subsystems
% -- blocked reactions
% -- blocked transport reactions
% -- reactions with E.C. number(s)
% -- metabolic reactions without E.C. number
% -- assigned CIDs
% -- number of mass-unbalanced reactions
% -- number of charge-unbalanced reactions
%
% > a file summarizing the number of reactions per subsystem 'subsystem_rxns.txt'
%
% > a file with the media components 'medium.txt'
%
% > a file summarizing the biomass reaction 'biomass_reaction.txt'
%
% > a file summarizing the compartmentation 'compartmentation.txt'

model.S = full(model.S);
% ~~~~~~~~~~~~~~~~ general stats ~~~~~~~~~~~~~~~~ %
disp('general stats')
%% transport reactions
t=ismember(model.rxns,findRxnsFromSubSystem(model,'Transport'));
n_transport = sum(t);
clear formulas

%% subsystems
subsystems = unique([model.subSystems{:}])';
n_subsystems = numel(subsystems);

%% blocked reactions
blocked = findBlockedReaction(model);
n_blocked = numel(blocked);
n_blocked_tr = numel(setdiff(blocked, model.rxns(t)));
clear blocked

%% reactions with E.C. association
rxns_no_ec = cellfun(@(x)isequal(x,'null')|isempty(x), model.rxnECNumbers);
n_ec = sum(~rxns_no_ec);

%% mass- and charge balancing
model.S = sparse(model.S);
[~, imBalancedMass, imBalancedCharge, imBalancedRxnBool, ~, ~] = ...
    checkMassChargeBalance(model);
imBalancedMass = ~cellfun(@isempty, imBalancedMass);
imBalancedCharge = imBalancedCharge~=0 | isnan(imBalancedCharge);
imBalancedMassCharge = imBalancedMass & imBalancedCharge;

fid = fopen(fullfile(out_dir, 'model_stats.txt'), 'w');
fprintf(fid, '%s\t%d\n',...
    'reactions',                numel(model.rxns),...
    'metabolites',              numel(model.mets),...
    'genes',                    numel(model.genes),...
    'gpr',                      sum(~(cellfun(@isempty, model.rules))),...
    'no_gpr_metabolic',         sum(~(cellfun(@isempty, model.rules(~t)))),...
    'transport',                n_transport,...
    'subsystems',               n_subsystems,...
    'blocked',                  n_blocked,...
    'blocked_transport',        n_blocked_tr,...
    'reactions_with_EC',        n_ec,...
    'no_EC_metabolic',          sum(rxns_no_ec&~t),...
    'metabolites_CID',          sum(~cellfun(@isempty, model.metPubChemID)),...
    'imbalanced_reactions',     sum(imBalancedRxnBool),...
    'mass_imbalanced',          sum(imBalancedMass),...
    'charge_imbalanced',        sum(imBalancedCharge),...
    'mass_charge_imbalanced',   sum(imBalancedMassCharge),...
    'irreversible_reactions',   sum(xor(model.lb,model.ub)));
fclose(fid);

model.S = full(model.S);

clear imBalancedMassCharge imBalancedCharge imBalancedMass imBalancedRxnBool ...
    t n_ec n_blocked n_blocked_tr n_subsystems n_transport rxns_no_ec
%% ~~~~~~~~~~~~~~~~ reactions per subsystem ~~~~~~~~~~~~~~~~ %
disp('subsystems')
rxns_per_subsystem = cellfun(@(x)sum(ismember([model.subSystems{:}], x)),...
    subsystems);

fid = fopen(fullfile(out_dir, 'subsystem_rxns.txt'), 'w');
for i=1:numel(subsystems)
    fprintf(fid, '%s\t%d\n', subsystems{i}, rxns_per_subsystem(i));
end
fclose(fid);
clear subsystems rxns_per_subsystem

%% ~~~~~~~~~~~~~~~~ medium ~~~~~~~~~~~~~~~~ %
disp('medium')
exchange_rxns = model.rxns(contains(model.rxnNotes, 'medium'));
medium = regexp(exchange_rxns, 'cpd\d+\[e0\]', 'match');
medium_ids = [medium{:}]';
medium_names = model.metNames(findMetIDs(model, medium_ids));
medium_names = regexprep(medium_names, '_plus', '+');
medium_names = regexprep(medium_names, '_', ' ');
medium_names = regexprep(medium_names, 'e0$', '');
medium_names = strtrim(medium_names);
medium_ids = strtok(medium_ids, '[');

fid = fopen(fullfile(out_dir, 'medium.txt'), 'w');
for i=1:numel(medium_ids)
    fprintf(fid, '%s\t%s\n', medium_ids{i}, medium_names{i});
end
fclose(fid);

clear exchange_rxns medium medium_ids medium_names

%% ~~~~~~~~~~~~~~~~ biomass ~~~~~~~~~~~~~~~~ %
disp('biomass reaction')
if sum(model.c) > 1
    warning('there are two objectives specified, chosing the first as the biomass reaction')
    idx = find(model.c);
    model.c(:) = 0;
    model.c(idx(1)) = 1;
end

biomass_ids = model.mets(model.S(:, model.c==1)~=0);
biomass_names = model.metNames(findMetIDs(model, biomass_ids));
biomass_names = regexprep(biomass_names, '_plus', '+');
biomass_names = regexprep(biomass_names, '_$', '');
biomass_names = regexprep(biomass_names, '_', '-');

biomass_coeff = model.S(findMetIDs(model, biomass_ids), model.c==1);
[biomass_coeff, order] = sort(biomass_coeff, 'ascend');
biomass_ids = biomass_ids(order);
biomass_names = biomass_names(order);

fid = fopen(fullfile(out_dir, 'biomass.txt'), 'w');
for i=1:numel(biomass_ids)
    fprintf(fid, '%s\t%s\t%d\n', biomass_ids{i}, biomass_names{i}, biomass_coeff(i));
end
fclose(fid);

%% ~~~~~~~~~~~~~~~~ compartmentation ~~~~~~~~~~~~~~~~ %
disp('compartmentation')
mets_per_comp = cellfun(@(x)numel(findMetFromCompartment(model, x)),...
    model.comps);
rxns_per_comp = cellfun(@(x)size(findRxnFromCompartment(model, x),1),...
    model.comps);
compnames = regexprep(model.compNames, '_', ' ');

fid = fopen(fullfile(out_dir, 'compartmentation.txt'), 'w');
fprintf(fid, 'Compartment\tReactions\tMetabolites\n');
for i=1:numel(compnames)
    fprintf(fid, '%s\t%d\t%d\n', compnames{i}, rxns_per_comp(i), mets_per_comp(i));
end
fclose(fid);

%% ~~~~~ metabolite BRITE classification ~~~~~ %%
% all
disp('metabolite BRITE classification')
metBRITE = model.metBRITE;
metBRITE(cellfun(@isempty,metBRITE)) = {'not assigned'};
metBRITE = strtok(metBRITE, ';');
uniqClasses = unique(metBRITE);
nPerClass = cellfun(@(x)sum(ismember(metBRITE,x)), uniqClasses);

fid = fopen(fullfile(out_dir, 'met-brite.txt'), 'w');
for i=1:numel(uniqClasses)
    fprintf(fid, '%s\t%d\n', uniqClasses{i}, nPerClass(i));
end
fclose(fid);

% lipids
metBRITE = model.metBRITE(contains(model.metBRITE, 'Lipids;'));
metBRITE = cellfun(@(x)strsplit(x,';'),metBRITE, 'un', 0);
for i=1:numel(metBRITE)
    br = strtrim(metBRITE{i});    
    if numel(br)>2
        metBRITE{i} = strjoin(br([1 2 3]), '; ');
    elseif numel(br)>1
        metBRITE{i} = strjoin(br([1 2]), '; ');
    else
        metBRITE{i} = strjoin(br(1), '; ');
    end
end
uniqClasses = unique(metBRITE);
nPerClass = cellfun(@(x)sum(ismember(metBRITE,x)), uniqClasses);

fid = fopen(fullfile(out_dir, 'met-brite-lipids.txt'), 'w');
for i=1:numel(uniqClasses)
    fprintf(fid, '%s\t%d\n', uniqClasses{i}, nPerClass(i));
end
fclose(fid);
% other
metBRITE = model.metBRITE(contains(model.metBRITE, 'Other'));
metBRITE = cellfun(@(x)strsplit(x,';'),metBRITE, 'un', 0);
for i=1:numel(metBRITE)
    br = strtrim(metBRITE{i});
    if numel(br)>1
        metBRITE{i} = strjoin(br([1 2]), '; ');
    else
        metBRITE{i} = strjoin(br(1), '; ');
    end
end
uniqClasses = unique(metBRITE);
nPerClass = cellfun(@(x)sum(ismember(metBRITE,x)), uniqClasses);

fid = fopen(fullfile(out_dir, 'met-brite-other.txt'), 'w');
for i=1:numel(uniqClasses)
    fprintf(fid, '%s\t%d\n', uniqClasses{i}, nPerClass(i));
end
fclose(fid);

end