% investigate the effect of myristate uptake in predicted growth using FBA
clear;clc
% the solver for COBRA FBA
changeCobraSolver('ibm_cplex', 'all',0);
% model file
modelFile = 'model/iRi1574.mat';
load(modelFile)
bioIdx = model.c==1;

%% limit the uptake of other carbon sources than palmitate and myristate
cUptRxnIdx = findRxnIDs(model,{'r1533_e0','r1534_e0','r1006_e0','r1002_e0'});
% find optimal growth value and set lower bound of biomass reaction
s = optimizeCbModel(model);
model.lb(bioIdx) = s.f;
% minimize the sum of uptake fluxes through additional carbon source uptake
% reactions
model.c(:) = 0;
model.c(cUptRxnIdx) = 1;
model.osenseStr = 'min';
s = optimizeCbModel(model);
model.lb(cUptRxnIdx) = s.x(cUptRxnIdx);

%% lipid uptake
palmitateUptIdx = findRxnIDs(model,'r1007_e0');
myristateUptIdx = findRxnIDs(model,'r1633_e0');

% determine range of palmitate uptake at optimal growth
range = [0 0];
s = optimizeCbModel(model);
disp('growth: 16:0')
disp(s.x(bioIdx))
fprintf('associated 16:0 uptake flux: %.4g\n', s.x(palmitateUptIdx));

% set lower bound for growth to optimal value
model.lb(bioIdx) = s.x(bioIdx);

% minimize 16:0 uptake
model.c(:) = 0;
model.c(palmitateUptIdx) = 1;
model.osenseStr = 'min';
s = optimizeCbModel(model);
range(1) = s.x(palmitateUptIdx);

% maximize 16:0 uptake
model.osenseStr = 'max';
s = optimizeCbModel(model);
range(2) = s.x(palmitateUptIdx);

disp('range of palmitate uptake at optimal growth:')
disp(range)

% reset objective
model.c(:) = 0;
model.c(bioIdx) = 1;
model.lb(bioIdx) = 0;

% growth with unlimited uptake of palmitate and myristate
disp('growth: 16:0 + 14:0')
myrUB = 1000;
model.ub(myristateUptIdx) = myrUB;
s = optimizeCbModel(model);
disp(s.f)

% growth with limited palmitate uptake and without myrsitate
disp('growth: 0.1*16:0')
model.ub(myristateUptIdx) = 0;
model.ub(palmitateUptIdx) = 0.1*range(1);
s = optimizeCbModel(model);
growthPalmRed = s.f;
disp(growthPalmRed)

% complementation of limited palmitate with myristate
disp('growth: 0.1*16:0 + 14:0')
model.ub(myristateUptIdx) = myrUB;
s = optimizeCbModel(model);
growthPalmRedMyr = s.f;
disp(growthPalmRedMyr)

% percent increase:
disp('percent increase')
disp(100*(1-growthPalmRed/growthPalmRedMyr))

disp('associated myristate uptake flux')
disp(s.x(myristateUptIdx))
