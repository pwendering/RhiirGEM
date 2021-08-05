% general validation of the model
clear
clc
% top-level directory
topDir = '';
% set solver for COBRA Toolbox
changeCobraSolver('ibm_cplex', 'all',0);

% load model
modelFile = [topDir, 'model/iRi1574.mat'];
load(modelFile)


%% remove unused genes and split reversible reactions
model = removeUnusedGenes(model);
model_irr = convertToIrreversible(model);

%{
%% kcats
% all kcats from BRENDA and SABIO-RK
maxKcatFile = fullfile(topDir,'data/kcats/kcat-reference-data.tsv');
% kcats associated to model reactions (will be created if it does not exist
% yet)
modelKcatsFile = fullfile(topDir, 'data/kcats/kcats_model.txt');
% associate maximum kcat per reaction
kcats = assign_kcats(model_irr, maxKcatFile, 'fungi');
% unknown kcats are assigned the median value of known values
kcats(kcats<=0) = median(kcats(kcats>0));
% convert to unit per hour
kcats = kcats * 3600; % [h^-1]
%}
load('kcats_irrev')
%% protein molecular weights
% UniProt results for R. irregularis
uniProtFile = fullfile(topDir,'data/uniprot.tab');
% calculate molecular weights for model proteins
[~, mw] = getAccountedProteinModel(model_irr, uniProtFile, []);
% convert from  g mol^-1 to g mmol^-1
mw = mw / 1000; % [g / mmol]

%% total protein content
C = 0.106; % g gDW^-1
% multiply C by the same factor that has been used to ensure numerical 
% stability in the other analyses
% C = 1.2*C;

%% lipid uptake
bioIdx = findRxnIDs(model_irr,'r0648_c0');
model_irr.osenseStr = 'max';
model_irr.c(:) = 0;
model_irr.c(bioIdx) = 1;
model_irr.lb(bioIdx) = 0;
palmitateUptIdx = findRxnIDs(model_irr,'r1007_e0');
model_irr.lb(palmitateUptIdx) = 0;
model_irr.ub(palmitateUptIdx) = 1000;
myristateUptIdx = findRxnIDs(model_irr,'r1633_e0');
model_irr.lb(palmitateUptIdx) = 0;
model_irr.ub(myristateUptIdx) = 0;
model_irr.ub(contains(model_irr.rxnNames,'ngam')) = 1000;
% determine range of palmitate uptake at optimal growth
palmitateUptRange = [0 0];

% optimize enzyme-constraint model
[~,problem] = enzymeFBA(model_irr,kcats,mw,C,bioIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params = struct;
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

fprintf('growth with only 16:0 as lipid in medium: %.3g\n', x(bioIdx))
fprintf('associated 16:0 uptake flux: %.4g\n', x(palmitateUptIdx));

tmp=cellfun(@(x)regexp(x,'\d+','match'),model_irr.rules,'un',0);
tmp=[tmp{:}];
nRxnProtPairs = numel(tmp);
startIdx = numel(model_irr.rxns)+numel(model_irr.genes) + nRxnProtPairs + 1;
y = solution.x(startIdx : startIdx+numel(model_irr.genes)-1);
notExpressedIdx = find(y==0);
for i=notExpressedIdx'
    rxnIdx = find(contains(model_irr.rules, ['x(' num2str(i) ')']));
    for j=rxnIdx
        if x(j)~=0
            fprintf('Non-zero flux for reaction %s.\n',model_irr.rxns(j));
        end
    end
end


% set lower bound for growth to optimal value
f = .99;
model_irr.lb(bioIdx) = f*x(bioIdx);

% minimize 16:0 uptake
model_irr.c(:) = 0;
model_irr.c(palmitateUptIdx) = 1;
model_irr.osenseStr = 'min';

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,palmitateUptIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

palmitateUptRange(1) = x(palmitateUptIdx);

% maximize 16:0 uptake
model_irr.osenseStr = 'max';

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,palmitateUptIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

palmitateUptRange(2) = x(palmitateUptIdx);

fprintf('range of 16:0 uptake at %d%% of optimal growth: %.3g - %.3g\n', 100*f,palmitateUptRange)

% reset objective
model_irr.c(:) = 0;
model_irr.c(bioIdx) = 1;
model_irr.lb(bioIdx) = 0;

% growth with unlimited uptake of palmitate and myristate
myrUB = 1000;
model_irr.ub(myristateUptIdx) = myrUB;

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,bioIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

fprintf('growth with unlimited 14:0 and 16:0 in medium: %.3g\n', x(bioIdx))


% growth with limited palmitate uptake and without myrsitate
model_irr.ub(myristateUptIdx) = 0;
model_irr.ub(palmitateUptIdx) = 0.1*palmitateUptRange(2);

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,bioIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

growthPalmRed = x(bioIdx);
fprintf('growth with 10%% of maximum 16:0 uptake: %.3g\n', growthPalmRed)
mediumUptNoMyr = [num2cell(x(cellfun(@(x)isequal(x,{'Medium'}),model_irr.subSystems))),...
    model_irr.rxnNames(cellfun(@(x)isequal(x,{'Medium'}),model_irr.subSystems))];

% complementation of limited palmitate with myristate
model_irr.ub(myristateUptIdx) = myrUB;

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,bioIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

growthPalmRedMyr = x(bioIdx);
fprintf('growth with 10%% of maximum 16:0 uptake + 14:0: %.3g\n', x(bioIdx))

mediumUptWithMyr = [num2cell(x(cellfun(@(x)isequal(x,{'Medium'}),model_irr.subSystems))),...
    model_irr.rxnNames(cellfun(@(x)isequal(x,{'Medium'}),model_irr.subSystems))];

disp([mediumUptNoMyr(:,1) mediumUptWithMyr])


% percent increase:
fprintf('==> percent increase: %2.2g%%\n',100*(1-growthPalmRed/growthPalmRedMyr))

fprintf('14:0 uptake flux value : %.3g\n', x(myristateUptIdx))
fprintf('16:0 uptake flux value : %.3g\n', x(palmitateUptIdx))

% clearvars -except modelFile bioIdx

%% Pi export at 99% of optimal growth
model_irr.ub(palmitateUptIdx) = 1000;
% piSinkIdx = findRxnIDs(model_irr,'r1570_e0');
piSinkIdx = findRxnIDs(model_irr,'r0019_e0_b'); % better use transport from c to e

piUptRange = [0 0];

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,bioIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

model_irr.lb(bioIdx) = 0.99*x(bioIdx);

% minimize Pi uptake
model_irr.c(:) = 0;
model_irr.c(piSinkIdx) = 1;
model_irr.osenseStr = 'min';

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,piSinkIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

piUptRange(1) = x(piSinkIdx);

% maximize Pi uptake
model_irr.osenseStr = 'max';

[~,problem] = enzymeFBA(model_irr,kcats,mw,C,piSinkIdx,true);
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = model_irr.osenseStr;
problem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
params.NumericFocus = 3;
params.IntegralityFocus = 1;
params.OutputFlag = 0;
solution = gurobi(problem, params);
x = solution.x(1:numel(model_irr.rxns));

piUptRange(2) = x(piSinkIdx);

disp('range of Pi export at optimal growth:')
disp(piUptRange)