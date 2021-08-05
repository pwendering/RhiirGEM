% summarize model characteristics for plotting
clear
clc

modelFile = 'model/iRi1574.mat';
load(modelFile)

out_dir = 'results/stats';

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
% > a file summarizing the compartmentation 'compartmentation.txt'

model.S = full(model.S);

%% ~~~~~~~~~~~~~~~~ general stats ~~~~~~~~~~~~~~~~ %
disp('general stats')

% transport reactions
t = cellfun(@(x)numel(unique(regexp(x,'\[..\]','match'))),...
    printRxnFormula(model,'rxnAbbrList',model.rxns,'printFlag',false))>1;
n_transport = sum(t);
clear formulas

% subsystems
subsystems = unique([model.subSystems{:}])';
n_subsystems = numel(subsystems);

% blocked reactions
% blocked = findBlockedReaction(model);
[minFlux,maxFlux] = fva(model);
blocked = model.rxns(minFlux==0&maxFlux==0);
n_blocked = numel(blocked);
n_blocked_tr = numel(setdiff(blocked, model.rxns(t)));
clear blocked

% reactions with E.C. association
rxns_no_ec = cellfun(@(x)isequal(x,'null')|isempty(x), model.rxnECNumbers);
n_ec = sum(~rxns_no_ec);

% mass- and charge balancing
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