% Simulate growth on single carbon sources
clear; clc
% top-level directory
topDir = '';
% set solver for COBRA Toolbox
changeCobraSolver('ibm_cplex', 'all',0);
% iRi1574 model path
rhiirModelFile = [topDir, 'model/iRi1574.mat'];
load(rhiirModelFile);

carbonSources = {
    'm0564[e0]' % D-glucose
    'm0588[e0]' % D-fructose
    'm0047[e0]' % raffinose
    'm1282[e0]' % melibiose
    'm0366[e0]' % D-xylose (Sugiura et al., PNAS Vol. 117, 2020)
    'm0028[e0]' % acetate (Pfeffer et al., Plant Physiol. Vol. 120, 1999)
    'm1331[e0]' % myristate (Sugiura et al., PNAS Vol. 117, 2020)
    'm0129[e0]' % glycine
    'm1335[e0]' % gycerol (Bago et al., Plant Physiol. Vol. 131, 2003)
    'm0239[c0]' % trehalose
};

% find minimum uptake of palmitate at optimal biomass
minFlux = fluxVariability(model,'rxnNameList', 'r1007_e0');
model.ub(findRxnIDs(model,'r1007_e0')) = minFlux;

% remove carbon sources except palmitate
toRemove = {
    'r1533_e0' % D-glucose uptake
    'r1534_e0' % D-fructose uptake
    'r1006_e0' % myo-inositol uptake
    'r1002_e0' % glycine uptake
    'r1633_e0' % myristate uptake
};

model = removeRxns(model,toRemove,'metFlag',0);

% Growth without adding any carbon source except palmitate
s_1 = optimizeCbModel(model);
fprintf('Growth without adding any carbon source except palmitate %.4g\n', s_1.f)

growth = nan(numel(carbonSources),1);

for i=1:numel(carbonSources)
    tmpModel = addReaction(model,['EX_' carbonSources{i}],...
        ['-> ' carbonSources{i}]);

    s_2 = optimizeCbModel(tmpModel);
    growth(i) = s_2.f;
end

% display results in the console
disp([model.metNames(findMetIDs(model,carbonSources)) num2cell(growth)])

% prepare barplot
% x-labels
carbonSources = [{'m0384[e0]'}; carbonSources];
xNames = categorical(strrep(model.metNames(findMetIDs(model,carbonSources)),'_','-'));
xNames = reordercats(xNames,strrep(model.metNames(findMetIDs(model,carbonSources)),'_','-'));

% y-data
growth = [s_1.f; growth];

% barplot
bar(xNames,...
    growth,'black')
set(gca, 'LineWidth', 1.5,  'TickLength', [0.003 0])
ylabel('predicted growth rate [h^{-1}]', 'FontSize', 12)
% save figure as PNG
print([topDir 'results/figures/Figure-S3.png'], '-painters', '-dpng')