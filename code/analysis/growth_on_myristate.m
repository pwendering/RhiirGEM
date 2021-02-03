% general validation of the model
clear
clc

changeCobraSolver('ibm_cplex', 'all',0);
modelFile = 'model/iRi1572.mat';
load(modelFile)

%% lipid uptake

bioIdx = model.c==1;
palmitateUptIdx = findRxnIDs(model,'r1007_e0');
myristateUptIdx = findRxnIDs(model,'r1633_e0');

% determine range of palmitate uptake at optimal growth
range = [0 0];
s = optimizeCbModel(model);
disp('growth: 16:0')
disp(s.x(bioIdx))
fprintf('associated 16:0 uptake flux: %.4g\n', -s.x(palmitateUptIdx));

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
disp(-range)

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
model.lb(palmitateUptIdx) = 0.1*range(2);
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

clearvars -except modelFile bioIdx

%% Pi export
load(modelFile)

piSinkIdx = findRxnIDs(model,'r1570_e0');
range = [0 0];

s = optimizeCbModel(model);
model.lb(bioIdx) = s.x(bioIdx);

% minimize Pi uptake
model.c(:) = 0;
model.c(piSinkIdx) = 1;
model.osenseStr = 'min';
s = optimizeCbModel(model);
range(1) = s.x(piSinkIdx);

% maximize Pi uptake
model.osenseStr = 'max';
s = optimizeCbModel(model);
range(2) = s.x(piSinkIdx);

disp('range of Pi export at optimal growth:')
disp(range)