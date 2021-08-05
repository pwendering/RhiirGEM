% Determine important carbon sources under different fungal structure.
% upper bounds for reactions are limited by v_max which is estimated from
% transcript abundances and k_cat values
clear; clc

changeCobraSolver('ibm_cplex', 'all', 0);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('|                       START                        |')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('')
try
    parpool(4);
catch
    disp('Parallel pool already available')
end

topDir = '';

% new sampling or only analysis of results
new = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~~~~~~~~~~~~~~~~~ carbon sources ~~~~~~~~~~~~~~~~~ %
carbonSources = {
    'm0564[e0]' % D-Glucose
    'm0588[e0]' % D-Fructose
    'm0047[e0]' % Raffinose
    'm1282[e0]' % Melibiose
    };

monosaccharideImportRxns = {'r0822_e0_f', 'r0824_e0_f', 'r1531_e0_f'};
concentration = 1000;

% ~~~~~~~~~~~~~~~~~~~~~~ model ~~~~~~~~~~~~~~~~~~~~~~ %
modelFile = fullfile(topDir,'model/iRi1574.mat');

% uptake reactions to be removed
oldUptake = {
    'r1533_e0' % D-Glucose
    'r1534_e0' % D-Fructose
    'r1006_e0' % myo-Inositol
    'r1002_e0' % Glycine
    'r1633_e0' % Myristate
    };

% ~~~~~~~~~~~~~~ transcriptomics data ~~~~~~~~~~~~~~~ %
transcriptomicDir = fullfile(topDir, 'data/transcriptomic-data');
transcriptomicFiles = {...
    fullfile(transcriptomicDir, 'transcriptomic_ERM.csv'),...
    fullfile(transcriptomicDir, 'transcriptomic_HYP.csv'),...
    fullfile(transcriptomicDir, 'transcriptomic_ARB.csv'),...
    };
nStructures = numel(transcriptomicFiles);
experiments = regexp(transcriptomicFiles,'(ARB)|(ERM)|(HYP)', 'match');
experiments = [experiments{:}];
clear transcriptomicDir

% ~~~~~~~~~~~~~~~~~~~~~~ kcats ~~~~~~~~~~~~~~~~~~~~~~ %
maxKcatFile = fullfile(topDir, 'data/kcats/kcat-reference-data.tsv');
modelKcatsFile = fullfile(topDir, 'data/kcats/kcats_model.txt');

% ~~~~~~~~~~~~~~~~~~ UniProt data ~~~~~~~~~~~~~~~~~~~ %
uniProtFile = fullfile(topDir, 'data/uniprot.tab');

% ~~~~~~~~~~~~~~~~ output directory ~~~~~~~~~~~~~~~~~ %
outDir = fullfile(topDir, 'results/fungal-structures/');

% ~~~~~~~~~~~~~~~~ experimental data ~~~~~~~~~~~~~~~~ %
% set total protein content (maximum over all conditions [Hildebrand et al.
C = 0.106;
C_sd = 3.5*1E-3; % sd for max protein content;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~~~~~~~~~~~~~~~~~~~~~~ model ~~~~~~~~~~~~~~~~~~~~~~ %
fprintf('\nPreparing the model\n')
% load the model file
load(modelFile)

% limit palmitate uptake to the minimum flux at optimal growth using
% standard FBA
palmitateIdx = findRxnIDs(model,'r1007_e0');
tmpModel = model;
% limit biomass reaction to optimal growth rate
s = optimizeCbModel(tmpModel);
tmpModel.lb(model.c==1) = (1-1e-6)*s.f;
% change objective to minimization of palmitate uptake
tmpModel.c(:) = 0;
tmpModel.c(palmitateIdx) = 1;
tmpModel.osenseStr = 'min';
s = optimizeCbModel(tmpModel);
% set palmitate uptake to the minimum possible flux at optimal growth
model.ub(palmitateIdx) = s.f;
clear tmpModel s palmitateIdx

% remove blocked reactions
% add additional carbon sources temporarily to reactions acting on it are not blocked
model = addReaction(model, 'tmpRxn', 'reactionFormula', ['-> ' strjoin(carbonSources,' + ')],...
    'printLevel', 0);
[minFlux, maxFlux] = fva(model,1-1e-6);
blockedReactions = model.rxns(minFlux==0&maxFlux==0);
% blockedReactions = findBlockedReaction(model);
model = removeRxns(model, blockedReactions, 'metFlag', false);
model = removeRxns(model, 'tmpRxn', 'metFlag', false);

% split reversible reactions
model = convertToIrreversible(model);

% remove genes which are not used in any reaction
model = removeUnusedGenes(model);
model = rmfield(model,'rxnGeneMat');

% ~~~~~~~~~~~~~~~~~~~~~~ kcats ~~~~~~~~~~~~~~~~~~~~~~ %
fprintf('\n> associating kcats to reactions\n')

% assign reaction-specific kcats from EC numbers associated to the
% reactions; read from modelKcatsFile if already done before
kcats = assign_kcats(model, maxKcatFile, 'fungi', modelKcatsFile);

% assign the median kcat of non-zero values to all reactions with unknown
% kcat; lower of equal to zero because some values are negative as obtained
% from BRENDA
kcats(kcats<=0) = median(kcats(kcats>0));

% multiply with 3600 to obtain [h^-1] as unit
kcats = kcats * 3600;

% remove previously uptake reactions for carbon sources
idxRemove = contains(model.rxns, oldUptake);
model = removeRxns(model, model.rxns(idxRemove));

% remove kcats of removed reactions
kcats(idxRemove) = [];
clear idxRemove

% ~~~~~~~~~~~~~~ add uptake reactions ~~~~~~~~~~~~~~~ %
% add uptake reactions for all carbon sources
fprintf('> adding uptake reactions\n')
uptakeRxns = repmat({''}, numel(carbonSources),1);
cSourceNames = repmat({''}, numel(carbonSources),1);
nRxns = numel(model.rxns);
for i=1:numel(carbonSources)
    cSourceNames{i} = model.metNames{findMetIDs(model, carbonSources(i))};
    uptakeRxns{i} = ['r',num2str(nRxns+1),'_e0'];
    
    model = addReaction(model, uptakeRxns{i},...
        'reactionName', [cSourceNames{i}, '_uptake'],...
        'reactionFormula', ['-> ', carbonSources{i}],...
        'upperBound', concentration);
    nRxns = numel(model.rxns);
    
end
clear metName oldUptake
model.subSystems(end-3:end) = {{'Medium'}};
uptakeIdx = findRxnIDs(model, uptakeRxns);

% limit the transport/influx reactions for the respective carbon
% sources depending of the concentrations (hexose transporter kinetics
% were taken from S. cerevisiae as not all values were available for R. irregularis)
% > use high affinity transport kinetics
% > use lowest Km and highest Vmax measured to avoid
%   overconstraining the model (+/- standard deviation,
%   respectively)
% sources: Reifenberger et al. 1997, European Journal of Biochemistry, 245: 324-333
%          Meijer et al. 1996, Biochimica et Biophysica Acta (BBA) - Bioenergetics, 1277(3), 209-216
%               Km (mM)              Vmax (nmol min^-1 mg^-1)
% glucose:        0.7                       173
% fructose:       8.7                       15.6 * (1000 / 60) = 260
% galactose:      0.8                       139
% v = (Vmax * [S]) / (Km + [S])
% concentrations are given in [mM]
vMaxFactor = 60 / 1000;

% find transport reaction indices for carbon sources or the respective
% monosaccharides, which can be imported (forward direction)
transRxnIdx = ismember(model.rxns, monosaccharideImportRxns);
if sum(transRxnIdx)~=numel(monosaccharideImportRxns)
    error('at least one import reaction ID was not found in the model')
end

% get maximum for galactose import
% Raffinose
MWRaf = 504.46; % [g/mol]
MWGal = 180.156; % [g/mol]
MWSuc = 342.2965; % [g/mol]
wv = (MWRaf/1000)*concentration; % [g/l]
ratioGal = (MWGal/(MWSuc+MWGal));
concGal = (1000*ratioGal*wv)/MWGal; % [mM]
uptakeGalRaffinose = (vMaxFactor * 139 * concGal) / (0.8 + concGal);
% Melibiose
MWMel = 342.297; % [g/mol]
MWGal = 180.156; % [g/mol]
wv = (MWMel/1000)*concentration; % [g/l]
concGal = (1000*0.5*wv)/MWGal; % [mM]
uptakeGalMelibiose = (vMaxFactor * 139 * concGal) ./ (0.8 + concGal);

clear concGal MWRaf MWGal MWSuc wv ratioGal concGal MWMel
uptakeFluxes = [
    (vMaxFactor * 173 * concentration) / (0.7 + concentration);... D-Glucose
    (vMaxFactor * 260 * concentration) / (8.7 + concentration);... D-Fructose
    max(uptakeGalRaffinose,uptakeGalMelibiose)... maximum for Gal from Raffinose and Melibiose
];

clear vMaxFactor uptakeGalRaffinose uptakeGalMelibiose

% assign upper bounds on import fluxes
model.ub(transRxnIdx) = uptakeFluxes;

% add zero kcats for the uptake reactions
kcats(end:end+numel(carbonSources)) = 0;

% ~~~~~~~~~~~~~~~ sampling iterations ~~~~~~~~~~~~~~~ %
n = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new
    fprintf('\n> starting with growth simluation using transcriptomic data sets\n\n')
    
    % initialize matrix for optimal growth rates for mean, mean+sd, mean-sd
    growthMat = zeros(3,nStructures);
    
    for i=1:nStructures
        
        
        disp('~~~~~~~~~~~~~~~~~~~~~~~~')
        fprintf('Current structure: %s\n', experiments{i})
        disp('~~~~~~~~~~~~~~~~~~~~~~~~')
        
        % fraction of total protein mass that is is accounted for in the model
        [f, MW] = getAccountedProteinModel(model, uniProtFile, transcriptomicFiles{i});
        fprintf('Accounted mass fraction: %.3f\n', f)
        f = 1;
        fprintf('f set to %.1f\n', f)
        
        % ~~~~~~~~~~~~~~~~~~~~~~ 1. FBA ~~~~~~~~~~~~~~~~~~~~~~ %
        disp('> Flux Balance Analysis')
        
        % +- standard deviation: 
        % rescale biomass coefficients to protein content
        C_tmp = C + C_sd;
        tmpModel = rescaleBiomassCoefficients(model,'m0995[c0]',C_tmp*f);
        growthPlusSd = simulateGrowthVmax(tmpModel, kcats, transcriptomicFiles{i},...
            C_tmp, f, MW);
        
        % rescale biomass coefficients to protein content
        C_tmp = C - C_sd;
        tmpModel = rescaleBiomassCoefficients(model,'m0995[c0]',C_tmp*f);
        growthMinusSd = simulateGrowthVmax(tmpModel, kcats, transcriptomicFiles{i},...
            C_tmp, f, MW);
        clear C_tmp tmpModel
        
        % mean protein content:
        
        % rescale biomass coefficients to protein content
        model = rescaleBiomassCoefficients(model,'m0995[c0]',C*f);

        [growth,v,ecModel] = simulateGrowthVmax(model, kcats, transcriptomicFiles{i},...
            C, f, MW);
        
        growthMat(:,i) = [growth; growthMinusSd; growthPlusSd];
        
        % ~~~~~~~~~~~~~~~~~~~~~~ 2. FVA ~~~~~~~~~~~~~~~~~~~~~~ %
        disp('> Flux Variability Analysis')        
        [minFlux, maxFlux] = fva(ecModel);
        
        % ~~~~~~~~~~~~~~~~~~~ 3. sampling ~~~~~~~~~~~~~~~~~~~~ %
        disp(['> Sampling ',num2str(n),' flux distributions'])
        fluxSamples = distSampling(ecModel,n,minFlux,maxFlux);

        % ~~~~~~~~~~~~~~~~ 4. write results ~~~~~~~~~~~~~~~~~~ %
        
        % write FBA and FVA results to file
        writetable(array2table([v,minFlux,maxFlux],'RowNames',model.rxns,...
            'VariableNames', {'fbaFlux','minFlux','maxFlux'}),...
            [outDir,'fva-',experiments{i},'.csv'],...
            'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t');
        
        % write sampling results to file
        writetable(array2table(fluxSamples,'RowNames',model.rxns,...
            'VariableNames',strcat('iter_',strtrim(cellstr(num2str((1:n)'))))),...
            [outDir,'samples-',experiments{i},'.csv'],...
            'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t');
        
        clear minFlux maxFlux growth fluxSamples ecModel v
        
        fprintf('\n')
    end; clear ecModel f MW C
    
    % write matrix with growth rates to file
    writetable(array2table(growthMat,'RowNames',{'av','av+sd','av-sd'},...
            'VariableNames',experiments),...
            [outDir,'growth-rates-structures.csv'],...
            'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fluxSamplingMat = zeros(nRxns,n*nStructures);
colIdx=1:n;

% read sampling results
for i=1:nStructures
    tmpResTab = readtable([outDir,'samples-',experiments{i},'.csv'],...
        'ReadVariableNames',true,'ReadRowNames',true,'Delimiter','\t');
    fluxSamplingMat(:,colIdx+(i-1)*n) = table2array(tmpResTab);
end
clear tmpResTab

% ~~~ find reactions with large effect size between structures ~~~ %
% calculate non-parametric estimator for common language A_w
A_w = zeros(nRxns,0.5*nStructures*(nStructures-1));

for rxnIdx=1:nRxns
    c=0;
    for i=1:nStructures-1
        for j=i+1:nStructures
            c=c+1;
            p = fluxSamplingMat(rxnIdx,(i-1)*n+1:i*n);
            q = fluxSamplingMat(rxnIdx,(j-1)*n+1:j*n);
            
            A_w(rxnIdx,c) = sum(arrayfun(@(x)sum(x>q) + .5*sum(x==q),p)) / n^2;
        end
    end
    if mod(i,500)==0
        fprintf('calculated A_w for %d reactions\n', i)
    end
end

% write effect sizes to file
writetable(cell2table([model.rxns, model.rxnNames,...
    vertcat(model.subSystems{:}), num2cell(A_w)]),...
    [outDir, 'effect-size-table.tsv'],'WriteVariableNames',false,...
    'FileType','text')


largeEffectIdx = A_w>0.8;
largeEffectRxns = repmat({'-'},max(sum(largeEffectIdx)),size(A_w,2));
for i=1:size(A_w,2)
%     tmpNames = unique(model.rxnNames(largeEffectIdx(:,i)));
    tmpNames = unique(model.rxnKEGGID(largeEffectIdx(:,i)));
    largeEffectRxns(1:numel(tmpNames),i)=tmpNames;
end

% write reaction names with large effect size to file
writetable(cell2table(largeEffectRxns,'VariableNames',{'ERM_vs_HYP',...
    'ERM_vs_ARB', 'HYP_vs_ARB'}),...
    [outDir, 'large-effect-rxns.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t');

% find reactions with large effect sizes wrt to fluxes between structures

fprintf('total ERM - IRM: %d\n', sum(~ismember(largeEffectRxns(:,1),'-')))
fprintf('total ERM - ARB: %d\n', sum(~ismember(largeEffectRxns(:,2),'-')))
fprintf('total IRM - ARB: %d\n', sum(~ismember(largeEffectRxns(:,3),'-')))

fprintf('unique ERM - IRM: %d\n', numel(setdiff(largeEffectRxns(:,1), largeEffectRxns(:,[2 3]))))
fprintf('unique ERM - ARB: %d\n', numel(setdiff(largeEffectRxns(:,2), largeEffectRxns(:,[1 3]))))
fprintf('unique IRM - ARB: %d\n', numel(setdiff(largeEffectRxns(:,3), largeEffectRxns(:,[1 2]))))

fprintf('intersect comp 1 -- 2 : %d\n', numel(intersect(largeEffectRxns(:,1),largeEffectRxns(:,2))))
fprintf('intersect comp 1 -- 3 : %d\n', numel(intersect(largeEffectRxns(:,1),largeEffectRxns(:,3))))
fprintf('intersect comp 2 -- 3 : %d\n', numel(intersect(largeEffectRxns(:,2),largeEffectRxns(:,3))))

fprintf('intersect all : %d\n', numel(intersect(largeEffectRxns(:,1),...
    intersect(largeEffectRxns(:,2),largeEffectRxns(:,3)))))

%% investigate large-effect reactions using different threshold for A_w
t = [0.6 0.7 0.8];
uniqSubSyst = unique([model.subSystems{:}],'stable');
comps = {'ERM_vs_HYP', 'ERM_vs_ARB', 'HYP_vs_ARB'};
plotDim = [numel(t) numel(comps)];
% define colors for color map (RBG as in figure 3)
pieColors = [[140,81,61];
            [111,196,95];
            [207,104,185];
            [131,208,180];
            [205,79,110];
            [199,201,81];
            [119,69,179];
            [83,109,63];
            [108,120,200];
            [207,90,50];
            [211,166,102];
            [140,164,199];
            [113,74,112]]/255;

close all
fig = figure;
colormap(fig, pieColors);
c=0;
dummyLabels = repmat({''},size(uniqSubSyst));

for i=1:numel(t)
    for j=1:numel(comps)
        c=c+1;
        % define position
        subplot(plotDim(1),plotDim(2),c)
        rxnIdx = A_w(:,j)>=t(i);
        if sum(rxnIdx)>0
            % find subsystems associated with current reactions
            subsystems = [model.subSystems{rxnIdx}];
            % determine the number of reactions per subsystem (complete
            % list)
            nPerSubSyst = cellfun(@(x)sum(ismember(subsystems,x)),uniqSubSyst);
            % plot pie chart with empty labels
            pie(nPerSubSyst,dummyLabels)
            % add number of reactions above threshold below the pie chart
            text(0.5,-.1,['n = ' num2str(sum(rxnIdx))], 'units', 'normalized',...
                'HorizontalAlignment', 'center', 'FontSize', 14);
        else
            % of no reactions are found, plot an invisible pie chart
            p = pie([0 1]);
            set(p,'Visible','off')
        end
        
        if j==1
            % in the beginning of each row, display the threshold
            text(-0.8,0.5,['A_w \geq ' num2str(t(i))],'units', 'normalized',...
                'FontSize', 16)
        end
        
        % alter the outer margins of the subplots
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        ti(:) = 0;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];

        % add the title (pairwise comparison) above each column of pies
        if i==1
            title(strrep(comps(j),'_',' '), 'FontSize', 16)
        end
    end
end
% add a legend
legend(uniqSubSyst, 'Box', 'off', 'Position',...
    [0.38 0.08 0.19 0.25], 'FontSize',8, 'NumColumns', 2)
% save a PNG
print('large_effect_reactions.png','-painters','-dpng')

%% investigate which reactions have very large or very small distances to their respective upper bounds
% fraction of protein accounted for in the model (set to 1 even though this is not the case)
f = 1;
fig = figure;
colormap(fig, pieColors);
dummyLabels = repmat({''},size(uniqSubSyst));

% loop over structures
for i=1:numel(transcriptomicFiles)
    
    % create an enzyme-constraint model to new upper bounds
    [~, MW] = getAccountedProteinModel(model, uniProtFile,[]);
    [~,~,ecModel] = simulateGrowthVmax(model, kcats, transcriptomicFiles{i},...
        C, f, MW);
    UB = ecModel.ub;
    
    % get average flux from sampling
    v = mean(fluxSamplingMat(:,(i-1)*n+1:i*n),2);
    
    UB(v==0) = NaN;
    v(v==0) = NaN;
    
    % pie chart for reactions that have a large distance to their upper
    % bound
    q90 = quantile(UB-v,0.9);
    subsystems = [model.subSystems{(UB-v)>=q90}];
    nPerSubSyst = cellfun(@(x)sum(ismember(subsystems,x)),uniqSubSyst);
    
    subplot(2,4,i)
    pie(nPerSubSyst,dummyLabels)
    title(experiments{i})
    if i==1
        text(-0.5,0.5,'>= q_{90}(UB-v)','units', 'normalized',...
                'FontSize', 14)
    end
    
    % pie chart for reactions that have a small distance to their upper
    % bound
    q10 = quantile(abs(UB-v),0.1);
    subsystems = [model.subSystems{abs(UB-v)<q10}];
    nPerSubSyst = cellfun(@(x)sum(ismember(subsystems,x)),uniqSubSyst);
    
    subplot(2,4,i+numel(transcriptomicFiles)+1)
    pie(nPerSubSyst,dummyLabels)
    
    % display considered fraction in the beginning of each row
    if i==1
        text(-0.5,0.5,'< q_{10}(UB-v)','units', 'normalized',...
            'FontSize', 14)
    end
    
    % calculate the overlap of these reaction with large-effect reactions
    largeEffectIdx = any(A_w>=0.6,2);
    fprintf('Average distance of large effect reaction to upper bounds: %.2g\n',...
        mean(abs(UB(largeEffectIdx)-v(largeEffectIdx))))
    fprintf('Average distance of all reactions to upper bounds: %.2g\n',...
        mean(abs(UB-v),'omitnan'))
end
% add a legend (subsystems)
legend(uniqSubSyst, 'Box', 'off', 'Position',...
    [0.78 0.41 0.18 0.22])
% save figure
print('rhiir_subsyst_rxns_large_diff_ub.png', '-painters', '-dpng')