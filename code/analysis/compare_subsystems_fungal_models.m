% Compare EC number overlab and fluxes through YeastGEM v8.3.4 subsystems accross
% published fungal models
clear
clc
% top-level directory
topDir = '';
% set solver for COBRA Toolbox
changeCobraSolver('ibm_cplex', 'all',0);
% iRi1572 model path
rhiirModelFile = [topDir, 'model' filesep 'iRi1574.mat'];
% path to EC correction file
ec_path = [topDir 'data' filesep 'corrected-EC-numbers.csv'];

fungalModelDir = [topDir 'results' filesep 'fungal_models'];

fungalModelFilenames = {
    'P_Chrysogenum_iAL1006_v1_00.xlsx' % works -- had to set objective, RAVEN style
    'P_pastoris_iRY1243.xlsx'   % failed -- file cannot be read
    'yeastGEM.xml'              % works
    'R_toruloides.xml'          % works
    'A_gossypii_iRL766.sbml'    % works
    'L_kluyveri_iPN730.xml'     % works but no EC numbers
    'M_alpina_iCY1106.xml'      % works
    'K_lactis_iOD907.xml'       % works
    'A_terreus_iJL1454_modified.xlsx'    % works -- had to re-format excel spreadsheet and remove 5 exchange reactions
    'Y_lipolytica_iYL_2_0_modified.xlsx' % failed -- added definitions for last 9 metabolites
    'N_crassa_iJDZ836.xml'      % works
};

modelNames = {
    '{\it P. chrysogenum} (iAL1006)'
    '{\it P. pastoris} (iRY1243)'
    '{\it S. cerevisiae} (YeastGEM v8.3.5)'
    '{\it R. toruloides} (rhto-GEM v1.3.0)'
    '{\it A. gossypii} (iRL766)'
    '{\it L. kluyveri} (iPN730)'
    '{\it M. alpina} (iCY1106)'
    '{\it K. lactis} (iOD907)'
    '{\it A. terreus} (iJL1454)'
    '{\it Y. lipolytica} (iYL 2.0)'
    '{\it N. crassa} (iJDZ836)'
    '{\it R. irregularis} (iRi1574)'
};
models = cell(numel(fungalModelFilenames)+1,1);

fprintf('Reading models from file...\n')

% read models from file
models{1} = importExcelModel(fullfile(fungalModelDir,fungalModelFilenames{1}));
models{1}.c(findRxnIDs(models{1},'bmOUT')) = 1;
models{1}.rxnECNumbers = models{1}.eccodes;
models{1} = rmfield(models{1},'eccodes');

models{2} = {};

models{3} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{3}));

models{4} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{4}));

models{5} = importModel(fullfile(fungalModelDir,fungalModelFilenames{5}));
models{5}.c(findRxnIDs(models{5},'NBIOMASS')) = 1; % N-limited biomass reaction
models{5}.rxnECNumbers = models{5}.eccodes;
models{5} = rmfield(models{5},'eccodes');

models{6} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{6})); % no EC

models{7} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{7}));

models{8} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{8}));
models{8}.rxnECNumbers = cellfun(@(x)strjoin(regexp(x,'\d+\.\d+\.\d+\.[0-9\-]+','match'),'|'),...
    models{8}.rxnNotes, 'un',0);

% A. terreus modifications:
% removed exchange reactions "D-Erythrose-4-phosphate exchange",
% "Galactitol exchange", "Inositol 1,3,4,5,6-pentakisphosphate",
% "6-methysalicylate acid exchange", "D-Ribose-5-phosphate exchange"
models{9} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{9}));
models{9}.c(findRxnIDs(models{9},'Biomass')) = 1;

models{10} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{10}));

models{11} = readCbModel(fullfile(fungalModelDir,fungalModelFilenames{11}));
tmpModel = importModel(fullfile(fungalModelDir,fungalModelFilenames{11}));
models{11}.rxnECNumbers = tmpModel.eccodes;

% load R. irregularis GEM
load(rhiirModelFile)
models{end} = model;

clear tmpModel fungalModelDir model topDir rhiirModelFile

% remove models that are not suited for further analysis
models([2 6]) = [];
modelNames([2 6]) = [];
fungalModelFilenames([2 6]) = [];

% check growth
growth = nan(numel(models),1);
for i=1:numel(models)
    s = optimizeCbModel(models{i});
    growth(i) = s.f;
end; clear s

% group reactions by subsystems from the yeast model:
% create a translation table for EC number to subsystem
ec_translation_tab = readtable(ec_path, 'ReadVariableNames', false);
yeastGEM = models{contains(fungalModelFilenames,'yeast')};
yeastSubSystEC = [yeastGEM.subSystems correctEC(yeastGEM.rxnECNumbers, ec_translation_tab)];
yeastSubSystEC(cellfun(@isempty,yeastSubSystEC(:,1)),1) = {''};

% match EC numbers to subsystems in YeastGEM
for i=1:numel(models)
    tmpModel = models{i};
    % create field with empty entries
    tmpModel.yeastSubSystems = repmat({''},numel(tmpModel.rxns),1);
    % correct EC numbers to current standard
    tmpModel.rxnECNumbers = correctEC(tmpModel.rxnECNumbers, ec_translation_tab);
    for j=1:numel(tmpModel.rxns)
        % find EC numbers associated with each reaction
        tmpEC = regexp(tmpModel.rxnECNumbers{j},'\d+\.\d+\.\d+\.[0-9\-]+', 'match');
        % remove those EC numbers that are not contained in the translation
        % table
        tmpEC(~cellfun(@(x)sum(contains(yeastSubSystEC(:,2),x)),tmpEC)) = [];
        % if EC numbers are found, get associated subsystems
        if ~isempty(tmpEC)
            tmpSubSyst = cellfun(@(x)yeastSubSystEC{contains(yeastSubSystEC(:,2),x),:},tmpEC, 'un',0);
            tmpModel.yeastSubSystems{j} = unique([tmpSubSyst{:}]);
        end
    end
    models{i} = tmpModel;
    fprintf('Mapped subsystems for model: %s\n',modelNames{i});
end; clear ec_translation_tab

% find models for which the mapping was successful
idxValid = find(cellfun(@(x)isfield(x,'yeastSubSystems'),models));
% get a list of unique yeast subsystems
uniqSubSyst = unique([yeastGEM.subSystems{:}])';
uniqSubSyst(cellfun(@isempty,uniqSubSyst)) = [];
% combine subsystem fields from all models
combSubSyst = cellfun(@(x)[x.yeastSubSystems{:}],models(idxValid),'un',0);
combSubSyst = [combSubSyst{:}];
% determine the total number of reactions matched to the subsystems accross
% all fungal models
nTotPerSubSystem = cellfun(@(x)sum(ismember(combSubSyst,x)),uniqSubSyst);
% define a set of subsystems for which flux should be compared between the models 
subSystForComp = {
    'sce01230  Biosynthesis of amino acids'
    'sce01200  Carbon metabolism'
    'sce01212  Fatty acid metabolism'
    'sce00230  Purine metabolism'
    'sce00240  Pyrimidine metabolism'
    'sce00010  Glycolysis'
    'sce00020  Citrate cycle (TCA cycle)'
    'sce00220  Arginine biosynthesis'
    };

% initialize matrices for pFBA, minimum and maximum flux values per
% subsystem and model
fluxPerSubsystem = zeros(5000,numel(subSystForComp)*numel(idxValid));
maxFluxPerSubsystem = zeros(5000,numel(subSystForComp)*numel(idxValid));
minFluxPerSubsystem = zeros(5000,numel(subSystForComp)*numel(idxValid));

fprintf('Determining flux ranges per subsystem...\n')

for i=1:numel(idxValid)
    tmpModel = models{idxValid(i)};
    tmpModel.subSystems = tmpModel.yeastSubSystems;
    
    % pFBA
    tmpModel = convertToIrreversible(tmpModel);
    % set maximum upper bound to 1000 for comparability
    tmpModel.ub(tmpModel.ub==max(tmpModel.ub)) = 1000;
    s_1 = optimizeCbModel(tmpModel, 'max', 'one');
    
    % set the minimum for growth to 90% of the optimal value
    tmpBioIdx = find(tmpModel.c);
    tmpModel.lb(tmpBioIdx) = .9*s_1.x(tmpBioIdx);
    
    % for each subsystem of interest, maximize/minimize the sum of fluxes
    % of the associated reactions:
    for j=1:numel(subSystForComp)
        tmpRxns = findRxnsFromSubSystem(tmpModel,subSystForComp{j});
        
        % maximize flux through individual subsystems at 90% of optimal growth
        tmpModel.c(:) = 0;
        tmpModel.c(findRxnIDs(tmpModel, tmpRxns)) = 1;
        s_2 = optimizeCbModel(tmpModel);
        
        % minimize flux through individual subsystems at 90% of optimal growth
        tmpModel.osenseStr = 'min';
        s_3 = optimizeCbModel(tmpModel);
        tmpModel.osenseStr = 'max';
        
        fluxPerSubsystem(1:numel(tmpRxns),i+((j-1)*numel(idxValid))) = ...
            s_1.x(findRxnIDs(tmpModel,tmpRxns));
        maxFluxPerSubsystem(1:numel(tmpRxns),i+((j-1)*numel(idxValid))) = ...
            s_2.x(findRxnIDs(tmpModel,tmpRxns));
        minFluxPerSubsystem(1:numel(tmpRxns),i+((j-1)*numel(idxValid))) = ...
            s_3.x(findRxnIDs(tmpModel,tmpRxns));
    end
    fprintf('\tDone with %s\n', modelNames{i});
end
% set values lower or equal to zero to NaN
fluxPerSubsystem(fluxPerSubsystem<=0) = NaN;
fluxPerSubsystem(all(isnan(fluxPerSubsystem),2),:) = [];
maxFluxPerSubsystem(maxFluxPerSubsystem<=0) = NaN;
maxFluxPerSubsystem(all(isnan(maxFluxPerSubsystem),2),:) = [];
minFluxPerSubsystem(minFluxPerSubsystem<=0) = NaN;
minFluxPerSubsystem(all(isnan(minFluxPerSubsystem),2),:) = [];
clear ec_path combSubSyst growth i j ia maxFlux minFlux s_1 s_2 s_3 tmpBioIdx ...
    tmpModel tmpRxns tmpSubSyst yeastSubSystEC tmpEC

% save workspace
save('rhiir_comp')
% load('rhiir_comp')

%% Compare minimum, maximum and pFBA flux sums per subsystem across fungal models
close all
figure

maxSumPerSubSyst = sum(maxFluxPerSubsystem,1, 'omitnan');
minSumPerSubSyst = sum(minFluxPerSubsystem,1, 'omitnan');
pFBASumPerSubSyt = sum(fluxPerSubsystem,1,'omitnan');

% print average coefficients of variation for all sums across the compared
% subsystems
cvMinSum = arrayfun(@(i)std(minSumPerSubSyst(((i-1)*numel(models)+1):i*numel(models)))/...
    mean(minSumPerSubSyst(((i-1)*numel(models)+1):i*numel(models))),1:numel(subSystForComp));
fprintf('Average coefficient of variation for minimum flux sums: %.2g\n',...
    mean(cvMinSum))

cvMaxSum = arrayfun(@(i)std(maxSumPerSubSyst(((i-1)*numel(models)+1):i*numel(models)))/...
    mean(maxSumPerSubSyst(((i-1)*numel(models)+1):i*numel(models))),1:numel(subSystForComp));
fprintf('Average coefficient of variation for maximum flux sums: %.2g\n',...
    mean(cvMaxSum))

cvpFBASum = arrayfun(@(i)std(pFBASumPerSubSyt(((i-1)*numel(models)+1):i*numel(models)))/...
    mean(pFBASumPerSubSyt(((i-1)*numel(models)+1):i*numel(models))),1:numel(subSystForComp));
fprintf('Average coefficient of variation for pFBA flux sums: %.2g\n',...
    mean(cvpFBASum))

% define x-positions for vertical lines
breakSize = 3;
xpos = arrayfun(@(i)(i-1)*(numel(idxValid)+breakSize)+1:(i-1)*(numel(idxValid)+breakSize)+numel(idxValid),...
    1:numel(subSystForComp),'un',0);
labelPos = cellfun(@mean,xpos);
% concatenate all X-positions
xpos = [xpos{:}];

% define colors
colors = [[125,130,178];[142,153,204];[217,219,241];[130,150,100];...
    [203,220,191];[236,176,80];[185,181,144];[188,133,145];[155,110,112];...
    [124,38,73]]/255;
lineColors = repmat(colors,numel(subSystForComp),1);

% plot lines
arrayfun(@(i)...
    line([xpos(i) xpos(i)], [minSumPerSubSyst(i) maxSumPerSubSyst(i)],...
    'Color',lineColors(i,:),'LineWidth',8),...
    1:numel(maxSumPerSubSyst))
hold on
% plot triangles/points for maximum, minimum and pFBA sums
scatter(xpos,maxSumPerSubSyst,30,lineColors-.1,'v','filled')
scatter(xpos,minSumPerSubSyst,30,lineColors-.1,'^','filled')
scatter(xpos,pFBASumPerSubSyt,30,lineColors-.2,'o','filled')
hold off
% x-labels and limits
xticks(labelPos)
xticklabels(regexprep(subSystForComp,'sce\d+',''))
xtickangle(-30)

xlim([0 max(xpos)+1])

% y-label
ylabel('Sum of fluxes [mmol gDW^{-1} h^{-1}]')

% add the names of fungal models legend
legend(regexprep(modelNames,'_','\\_'),'NumColumns', 5,...
    'FontSize',12, 'Position',[0.15 0.94 0.75 0.05],'Box','off',...
    'FontWeight', 'bold', 'TextColor', [0 0 0]);
box on
% change graphical parameters and transform Y-axis to log-scale
set(gca,'FontSize',14, 'LineWidth',1.2, 'TickLength', [.005, .01],...
    'YScale', 'log')
set(gcf, 'units','normalized','OuterPosition',[-1.3 -.4 1.3 1.3]);
% save figure as PNG
print([topDir 'results/figures/Figure-S2.png'], '-painters', '-dpng');

%% EC overlap per subsystem between iRi1572 model and other fungal models
close all
% chose the subsystems with most associated reactions for this comparison
% (n>50% quantile)
compSubsyst = uniqSubSyst(nTotPerSubSystem>=quantile(nTotPerSubSystem,.5));
pwOverlap = nan(numel(compSubsyst,numel(idxValid)-1));
iRi1574 = models{end};
iRi1574.subSystems = iRi1574.yeastSubSystems;

% compare subsystems using Jaccard index
for i=1:numel(compSubsyst)
    % EC numbers in iRi1574 model associated with the current subsystem
    ec1 = iRi1574.rxnECNumbers(findRxnIDs(iRi1574,findRxnsFromSubSystem(iRi1574,compSubsyst(i))));
    ec1 = regexp(ec1,'\d+\.\d+\.\d+\.[0-9\-]+', 'match');
    ec1 = unique([ec1{:}]);
    
    for j=1:numel(models)-1 % -1 since iRi1574 is the last model
        tmpModel = models{j};
        % EC numbers associated with the current subsystem in the compared
        % model
        tmpModel.subSystems = tmpModel.yeastSubSystems;
        ec2 = tmpModel.rxnECNumbers(findRxnIDs(tmpModel,findRxnsFromSubSystem(tmpModel,compSubsyst(i))));
        ec2 = regexp(ec2,'\d+\.\d+\.\d+\.[0-9\-]+', 'match');
        ec2 = unique([ec2{:}]);
        
        % if both models have EC numbers associated with the subsystem,
        % calculate the Jaccard Index
        if ~isempty(ec1) && ~isempty(ec2)
            pwOverlap(i,j) = numel(intersect(ec1,ec2)) / numel(union(ec1,ec2));
        end
    end
    
end

% prepare the figure
figure
% sort subsystems by overlap and plot heatmap
[~,ia] = sort(sum(pwOverlap,2),'descend');
heatmap(regexprep(modelNames(1:end-1),'_','\\_'),...
    regexprep(compSubsyst(ia),'^sce\d+',''),pwOverlap(ia,:),...
    'Colormap', pink,'FontSize',8, 'FontName', 'Arial')
% add the annotation for the colorbar
annotation(gcf,'textarrow',...
    [.96 .96],[0.33 0.33],...
    'String','Jaccard Index with iRi1574 E.C. numbers',...
    'HeadStyle','none','LineStyle','none',...
    'FontSize',12,'FontName', 'Arial',...
    'TextRotation',-90);
% set plot dimensions
set(gca,'Position', [.4 .17 .45 .78])
set(gcf, 'units','normalized','OuterPosition',[0 0 .5 1]);
% save figure as PNG
print([topDir 'results/figures/Figure-S1.png'],'-painters','-dpng')

%% Percentage of reactions matched to subsystems per fungal model
percentRxnsInSubsyst = zeros(numel(models),1);
for i=1:numel(models)
    percentRxnsInSubsyst(i) = sum(~cellfun(@isempty,models{i}.yeastSubSystems))/numel(models{i}.rxns);
end
