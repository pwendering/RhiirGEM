% plot distribution of kcat values for the iRi1574 model and YeastGEM
% v8.3.4
clear
clc
% top-level directory
topDir = '';
changeCobraSolver('ibm_cplex', 'all',0);
% kcats per subsystem
modelFile = [topDir, 'model/iRi1574.mat'];
load(modelFile)

%% remove unused genes and split reversible reactions
model = removeUnusedGenes(model);
model_irr = convertToIrreversible(model);

maxKcatFile = fullfile(topDir,'data/kcats/kcat-reference-data.tsv');

% associate maximum kcat per reaction
kcats = assign_kcats(model_irr, maxKcatFile, 'fungi');

disp(['Median kcat: ' num2str(median(kcats(kcats>0)))])

%% Plot kcat values per subsystem
uniqSubSyst = unique([model.subSystems{:}]);
kcatsPerSubSyst = zeros(numel(kcats),numel(uniqSubSyst));
for i=1:numel(uniqSubSyst)
    tmpKcat = kcats(cellfun(@(x)isequal(x,uniqSubSyst(i)),model.subSystems));
    kcatsPerSubSyst(1:numel(tmpKcat),i) = tmpKcat;
end
% log-transform kcats
kcatsPerSubSyst_log = log10(kcatsPerSubSyst);
kcatsPerSubSyst_log(isinf(kcatsPerSubSyst_log)) = NaN;
uniqSubSyst(all(isnan(kcatsPerSubSyst_log),1)) = [];
kcatsPerSubSyst_log(:,all(isnan(kcatsPerSubSyst_log),1)) = [];

% prepare figure
figure
boxplot(kcatsPerSubSyst_log,'Labels',...
    strcat(uniqSubSyst', strcat(' (n =',{' '}, num2str(sum(~isnan(kcatsPerSubSyst_log))'),')')),...
    'PlotStyle', 'compact','Symbol', '.', 'Colors', 'k', 'OutlierSize', 8,...
    'Orientation', 'horizontal','LabelOrientation','horizontal')
xlabel('log_{10} k_{cat} [s^{-1}]')
box on
set(gca,'LineWidth',1.5,  'FontSize', 16)
set(gcf, 'units','normalized','OuterPosition',[0 0 .7 1]);
% save figure as PNG
print([topDir 'results/figures/kcats_per_subsystem.png'], '-painters', '-dpng')

% Compare kcats to current yeast model
yeastGEM = readCbModel;
% run GECKO, extract matched kcats and write them into a table
% ecYeastGEM = enhanceGEM(yeastGEM, 'COBRA');

% read the table with extracted kcat values
gecko_kcat_tab = readtable([topDir 'data/kcats/ecYeastGEM_kcats.tab'],...
    'FileType', 'text', 'Delimiter', '\t',...
    'ReadVariableNames', true, 'Format','%s\t%s\t%s');

% run kcat matching function that was used for the iRi1574 model:
% compartment information must be removed from metabolite names
yeastGEM.metNames = strtrim(strtok(yeastGEM.metNames,'['));
yeastGEM_irr = convertToIrreversible(yeastGEM);
kcats_yeast = assign_kcats(yeastGEM_irr,maxKcatFile,'Saccharomyces cerevisiae');

%% plot kcats per EC between yeast and R. irregularis model
% find intersection of EC numbers
ecRhiir = regexp(model.rxnECNumbers,'\d+\.\d+\.\d+\.\d+', 'match');
ecRhiir = [ecRhiir{:}];
ecYeast = regexp(yeastGEM.rxnECNumbers,'\d+\.\d+\.\d+\.\d+', 'match');
ecYeast = [ecYeast{:}];
ecIS = intersect(ecRhiir, ecYeast);

% find all kcats that are assigned to EC numbers in both models
maxKcatMatchRhiir = cellfun(@(x)max(kcats(contains(model.rxnECNumbers,x))),ecIS);
maxKcatMatchYeast_simple = cellfun(@(x)max(kcats_yeast(contains(yeastGEM.rxnECNumbers,x))),ecIS);

maxKcatMatchYeast_gecko = nan(size(ecIS));
for i=1:numel(ecIS)
    ecIdx = find(contains(gecko_kcat_tab.EC_Numbers, ecIS{i}));
    if ~isempty(ecIdx)
        ecKcats = nan(size(ecIdx));
        for j=1:numel(ecIdx)
            tmpECIdx = contains(strsplit(gecko_kcat_tab.EC_Numbers{ecIdx(j)}, ';'),ecIS{i});
            tmpKcats = strsplit(gecko_kcat_tab.kcats{ecIdx(j)},';');
            ecKcats(j) = max(str2double(tmpKcats(tmpECIdx)));
        end
        
        maxKcatMatchYeast_gecko(i) = max(ecKcats);
    end
end

% scatter plot
close all
figure
plot(maxKcatMatchRhiir, maxKcatMatchYeast_simple, 'k.', 'MarkerSize', 10)
hold on
plot(maxKcatMatchRhiir, maxKcatMatchYeast_gecko, 'kv', 'MarkerSize', 5)

xlabel('Maximum k_{cat} per E.C. number (iRi1574) [s^{-1}]')
ylabel('Maximum k_{cat} per E.C. number (YeastGEM v8.3.4) [s^{-1}]')

legend({['maximum per reaction (\rho_P = ' num2str(corr(maxKcatMatchYeast_simple',...
    maxKcatMatchRhiir'),'%.2g') ')'], ['GECKO (\rho_P = ' num2str(corr(maxKcatMatchYeast_gecko',...
    maxKcatMatchRhiir', 'rows', 'complete'),'%.2g') ')']},'Box', 'off',...
    'Location', 'northwest')

box on
set(gca,'yscale','log', 'xscale', 'log', 'FontSize', 12, 'LineWidth', 1.3)
set(gcf, 'OuterPosition', [353.6667 87.0000 436.6666 612.0000])

print([topDir 'results/figures/kcats_rhiir_vs_yeast_max_gecko.png'], '-painters', '-dpng')


