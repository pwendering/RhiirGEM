% Perform coupling analysis of fungal models
clear
clc
% top-level directory
topDir = '';
% set solver for COBRA Toolbox
changeCobraSolver('ibm_cplex', 'all',0);
% iRi1572 model path
rhiirModelFile = [topDir, 'model' filesep 'iRi1574.mat'];

fungalModelDir = [topDir 'results' filesep 'fungal-models'];

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

%% Read models form file
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

%% Assess coupling to biomass reaction
alpha = 90;
rangeTol = 1e-3;


couplingTypes = {'hard-coupled','soft-coupled','partially-coupled', 'uncoupled', 'blocked'};
couplingMat = nan(numel(couplingTypes),numel(models));

for i=1:numel(alpha)
    for j=1:numel(models)
        disp(regexp(modelNames{j},'\w\. \w+','match'))
        tmpOpt = optimizeCbModel(models{j}).f;
        
        try
            [minFlux,maxFlux] = fva(models{j},alpha(i)/100);
            % [minFlux,maxFlux] = fluxVariability(models{j},alpha(i));
            ranges = maxFlux - minFlux;
            range_bio = ranges(models{j}.c==1);

            couplingMat(1,j) = sum(...
                abs(minFlux)>rangeTol & ...
                abs(maxFlux)>rangeTol & ...
                ranges >= (range_bio - rangeTol) & ...
                ranges <= (range_bio + rangeTol));
            couplingMat(2,j) = sum(...
                abs(minFlux)>rangeTol & ...
                abs(maxFlux)>rangeTol & ...
                ranges < (range_bio - rangeTol));
            couplingMat(3,j) = sum(...
                ranges > (range_bio + rangeTol) & ...
                ranges < (2000-rangeTol));
            couplingMat(4,j) = sum(...
                ranges >= 2000-rangeTol);
            couplingMat(5,j) = sum(...
                minFlux <= rangeTol & ...
                maxFlux <= rangeTol);
        catch ME
            disp(ME.message)
        end
        
    end
end

couplingTypes(end) = [];
couplingMat_wo_blocked = couplingMat(1:end-1,:);

couplingMat_scaled = couplingMat_wo_blocked * diag(1./sum(couplingMat_wo_blocked));

colors = [[142,153,204];[130,150,100];[236,176,80];[217,219,241];[124,38,73]]/255;

labels = categorical(modelNames);
labels = reordercats(labels,modelNames);

b = bar(labels,100*couplingMat_scaled','stacked', 'Horizontal','on');

for i=1:numel(couplingTypes)
    b(i).FaceColor = colors(i,:);
end

legend(couplingTypes,'Location','northoutside','NumColumns',numel(couplingTypes),...
    'FontSize',14, 'Box', 'off')
xlabel('Percentage of flux-carrying reactions')
xlim([0 100])
xticks([])
box off
set(gca,'FontSize',14)