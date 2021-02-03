% estimate enzymes usage on different carbon sources
% initCobraToolbox(false);
changeCobraSolver('ibm_cplex', 'all');

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('|                       START                        |')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('')
try
    parpool(3);
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

% concentrations used by Hildebrandt et al. 2006, FEMS: 10, 100, and 1000 mM
concentrations = [10 100 1000];
colNames = {'c1','c2','c3'}; % for writing result files

% ~~~~~~~~~~~~~~~~~~~~~~ model ~~~~~~~~~~~~~~~~~~~~~~ %
modelFile = fullfile(topDir,'Model_objects/max_model_final/iRi1572.mat');

% uptake reactions to be removed
oldUptake = {
    'r1533_e0' % D-Glucose
    'r1534_e0' % D-Fructose
    'r1006_e0' % myo-Inositol
    'r1002_e0' % Glycine
    'r1633_e0' % Myristate
    };

% metabolite ID of protein biomass component
proteinID='m0995[c0]';

% ~~~~~~~~~~~~~~~~~~~~~~ kcats ~~~~~~~~~~~~~~~~~~~~~~ %
maxKcatFile = fullfile(topDir,'data/kcats/kcat-reference-data.tsv');

modelKcatsFile = fullfile(topDir,'data/kcats/kcats_model.txt');

% ~~~~~~~~~~~~~~~~~~ UniProt data ~~~~~~~~~~~~~~~~~~~ %
uniProtFile = fullfile(topDir,'data/uniprot.tab');

% ~~~~~~~~~~~~~~~~ output directory ~~~~~~~~~~~~~~~~~ %
outDir = fullfile(topDir,'results/carbon-sources');

% ~~~~~~~~~~~~~~~~ experimental data ~~~~~~~~~~~~~~~~ %
% protein mass per gram dry weight [g/gDW] (Hildebrandt et al. 2006, FEMS)
proteinContent = [...
    91 95 96;               % Glucose
    88 98 89;               % Fructose
    102 103 106;            % Raffinose
    83 86 75] * 1E-3;       % Melibiose

proteinContent = 1.2*proteinContent;

% mass fraction of protein accounted for in the model (as in GECKO,
% Sanchez et al. 2017)
f = 1;
proteinContent=f*proteinContent;
clear f pContentSdev

% ~~~~~~~~~~~~~~~ sampling iterations ~~~~~~~~~~~~~~~ %
n = 1000;

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
tmpModel.lb(model.c==1) = s.f;
% change objective to minimization of palmitate uptake
tmpModel.c(:) = 0;
tmpModel.c(palmitateIdx) = 1;
tmpModel.osenseStr = 'min';
s = optimizeCbModel(tmpModel);
% set palmitate uptake to the minimum possible flux at optimal growth
model.lb(palmitateIdx) = s.f;
clear tmpModel s palmitateIdx

% remove blocked reactions
% add additional carbon sources temporarily to reactions acting on it are not blocked
model = addReaction(model, 'tmpRxn', 'reactionFormula', ['-> ' strjoin(carbonSources,' + ')],...
    'printLevel', 0);
blockedReactions = findBlockedReaction(model);
model = removeRxns(model, blockedReactions, 'metFlag', false);
model = removeRxns(model, 'tmpRxn', 'metFlag', false);

% split reversible reactions
model = convertToIrreversible(model);

% remove genes which are not used by any reaction
model = removeUnusedGenes(model);
model = rmfield(model,'rxnGeneMat');
nGenes = numel(model.genes);

% metabolite names for carbon sources
cSourceNames = model.metNames(findMetIDs(model, carbonSources));

% ~~~~~~~~~~~~~~~~~~~~~~ kcats ~~~~~~~~~~~~~~~~~~~~~~ %
fprintf('\nAssociating kcats to reactions\n\n')

% assign reaction-specific kcats from EC numbers associated to the
% reactions; read from modelKcatsFile if already done before
kcats = assign_kcats(model, maxKcatFile, 'fungi', modelKcatsFile);

% assign the median kcat of non-zero values to all reactions with unknown
% kcat (lower or equal to zero because some values are negative as obtained
% from BRENDA)
kcats(kcats<=0) = median(kcats(kcats>0));
kcats = kcats * 3600; % [h^-1]

% add a zero for the added exchange reaction
kcats(end+1) = 0;

% remove previously uptake reactions for carbon sources
idxRemove = contains(model.rxns, oldUptake);
model = removeRxns(model, model.rxns(idxRemove));
nRxns = numel(model.rxns)+1; % + uptake reaction, which will be added in each iteration

% remove kcats of removed reactions
kcats(idxRemove) = [];

% ~~~~~~~~ molecular weight for every enzyme ~~~~~~~~ %
[~, mw] = getAccountedProteinModel(model, uniProtFile, []);
mw = mw / 1000; % [g / mmol]

clear uniProtFile oldUptake maxKcatFile idxRemove

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over carbon sources and uptake fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new
    uptakeFluxes = zeros(numel(carbonSources),numel(concentrations));
    growthMatEnzFBA = zeros(numel(carbonSources), numel(concentrations));
    growthMatFBA = zeros(numel(carbonSources), numel(concentrations));
    
    for i=1:numel(carbonSources)
        fprintf('\n########### %s ###########\n', cSourceNames{i})
        
        % add carbon uptake reaction with default upper bound (1000)
        uptakeID = ['r',num2str(size(model.S,2)+1),'_e0'];
        tmpModel = addReaction(model, uptakeID,...
            'reactionName', [cSourceNames{i}, '_uptake'],...
            'reactionFormula', ['-> ', carbonSources{i}],...
            'upperBound', 1000,'printLevel',0);
        
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
        
        if isequal(cSourceNames{i},'d_glucose')
            uptakeFluxes = (vMaxFactor * 173 * concentrations) ./ (0.7 + concentrations);
            limitRxnIdx = find(ismember(model.rxnNames,'d_glucose_transport_e0'));
            
        elseif isequal(cSourceNames{i},'d_fructose')
            uptakeFluxes = (vMaxFactor * 260 * concentrations) ./ (8.7 + concentrations);
            limitRxnIdx = find(ismember(model.rxnNames,'d_fructose_transport_e0'));
            
        elseif isequal(cSourceNames{i},'raffinose')
            
            % Here, the trisaccharide is broken down to sucrose and
            % galactose. As sucrose cannot be further broken down by
            % R.irregularis, we only limit the influx of galactose based on
            % the splitting of mass.
            MWRaf = 504.46; % [g/mol]
            MWGal = 180.156; % [g/mol]
            MWSuc = 342.2965; % [g/mol]
            wv = (MWRaf/1000)*concentrations; % [g/l]
            ratioGal = (MWGal/(MWSuc+MWGal));
            concGal = (1000*ratioGal*wv)/MWGal; % [mM]
            % concSuc=(1000*(1-ratioGal)*wv)/MWSuc; % [mM]
            
            uptakeFluxes = (vMaxFactor * 139 * concGal) ./ (0.8 + concGal);
            limitRxnIdx = find(ismember(model.rxnNames,'d_galactose_transport_e0'));
            clear MWRaf MWGal MWSuc wv ratioGal concGal concSuc
            
        elseif isequal(cSourceNames{i},'melibiose')
            MWMel = 342.297; % [g/mol]
            MWGal = 180.156; % [g/mol]
            MWGlc = 180.156; % [g/mol]
            wv = (MWMel/1000)*concentrations; % [g/l]
            concGal = (1000*0.5*wv)/MWGal; % [mM]
            concGlc = (1000*0.5*wv)/MWGlc; % [mM]
            
            uptakeFluxes = [(vMaxFactor * 139 * concGal) ./ (0.8 + concGal);
                (vMaxFactor * 173 * concGlc) ./ (0.7 + concGlc)];
            limitRxnIdx = cell2mat(cellfun(@(x)find(ismember(model.rxnNames,x)),...
                {'d_galactose_transport_e0','d_glucose_transport_e0'},'un',0));
            clear MWMel MWGal MWGlc wv concGal concGlc
        end
        limitRxnIdx = limitRxnIdx(1,:);
        clear vMaxFactor
        
        for j=1:numel(concentrations)
            fprintf('\n### %d mM ###\n', concentrations(j))
            
            % set upper bound for uptake reaction
            tmpModel.ub(limitRxnIdx) = uptakeFluxes(:,j);
            
            % set protein content
            currentC = proteinContent(i,j); % [g/gDW]
            
            % rescale biomass coefficients according to protein content
            finalModel = rescaleBiomassCoefficients(tmpModel,proteinID,currentC);
            
            % ~~~~~~~~~~~~~~~~~~~~~~ 1. FBA ~~~~~~~~~~~~~~~~~~~~~~ %
            
            % run FBA without enzyme constraints
            s = optimizeCbModel(finalModel);
            
            growthMatFBA(i,j) = s.x(finalModel.c==1);
            
            % run FBA with enzyme constraints
            s = enzymeFBA(finalModel, kcats, mw, currentC, [], true);
            
            % store predicted growth
            growthMatEnzFBA(i,j) = s.x(finalModel.c==1);
            
            % % ~~~~~~~~~~~~~~ 2. variability analysis ~~~~~~~~~~~~~ %
            fprintf('\n> Variablity analysis for enzyme concentrations\n')
            
            [minConc,maxConc] = enzymeFVA(finalModel, kcats, mw, currentC);
            % ~~~~~~~~~~~~~~~~~~~ 3. sampling ~~~~~~~~~~~~~~~~~~~~ %
            fprintf('\n> Sampling %d enzyme concentrations at optimal growth\n', n)
            
            [eConcMat, rxnFluxMat] = enzymeSampling(finalModel, minConc, maxConc,...
                kcats, mw, currentC, n);
            
            clear finalModel
            % ~~~~~~~~~~~~ 4. write results to file ~~~~~~~~~~~~~~ %
            fprintf('\n> Writing results to file')
            
            % table of concentrations
            resTable = [
                cell2table(model.genes,...
                'VariableNames', {'GeneID'}),...
                array2table(s.x(nRxns+1:nRxns+nGenes),...
                'VariableNames', {'fbaConc'}),...
                array2table(minConc,...
                'VariableNames', {'minConc'}),...
                array2table(maxConc,...
                'VariableNames', {'maxConc'}),...
                array2table(var(eConcMat,1,2),...
                'VariableNames', {'variance'}),...
                array2table(eConcMat,...
                'VariableNames', strcat('iter_',strtrim(cellstr(num2str((1:n)')))))];
            
            % write concentration table
            writetable(resTable, [outDir, cSourceNames{i}, '_', colNames{j},...
                '_concentration_sampling.csv'],...
                'writeVariableNames', true, 'Delimiter', '\t')
            
            % write table of fluxes
            writetable(array2table([s.x(1:nRxns),rxnFluxMat],'RowNames',tmpModel.rxns,...
                'VariableNames',['fbaFlux'; strcat('iter_',strtrim(cellstr(num2str((1:n)'))))]),...
                [outDir, cSourceNames{i}, '_', colNames{j},...
                '_flux_sampling.csv'],...
                'WriteVariableNames',true,'WriteRowNames',true,'Delimiter','\t');
            
            fprintf('\n-----------------------------\n')
        end
        clear minConc maxConc eConcMat rxnFluxMat resTable uptakeID
    end
    
    % write table of optimal growth at the different conditions
    
    % FBA
    writetable(array2table(growthMatFBA, 'VariableNames', colNames,...
        'RowNames', cSourceNames),[outDir,'growth-fba.csv'],...
        'writeVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t')
    
    % FBA with enzyme constraints
    writetable(array2table(growthMatEnzFBA, 'VariableNames', colNames,...
        'RowNames', cSourceNames),[outDir,'growth.csv'],...
        'writeVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t')
    
end

clear rxnFluxMat s cMax currentC kcats ProtPerDW relProtDW ...
    excIdx growthMat uptakeFluxes f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCond = numel(carbonSources)*numel(concentrations);
concSamplingMat = zeros(nGenes,n*nCond);
fluxSamplingMat = zeros(nRxns,n*nCond);
fbaFluxes = zeros(nRxns,nCond);

% add dummy reaction for the respective uptake reaction per condition
model = addReaction(model,'uptake', 'reactionFormula','-> carbonSource',...
    'subSystem', 'Medium');

% conter of conditions
c=0;
% column indices
colIdx=1:n;
for i=1:numel(carbonSources)
    
    for j=1:numel(concentrations)
        
        % read the results of protein abundance sampling
        tmpResTab = readtable([outDir, cSourceNames{i}, '_',...
            colNames{j}, '_concentration_sampling.csv'],...
            'ReadVariableNames',true,'ReadRowNames',false);
        % select only the sampling columns
        tmpResTab = table2array(tmpResTab(:,6:end));
        concSamplingMat(:,colIdx+(n*c)) = tmpResTab;
        
        % read the results from reaction flux sampling
        tmpResTab = readtable([outDir, cSourceNames{i}, '_', colNames{j},...
            '_flux_sampling.csv'],'ReadVariableNames',true,...
            'ReadRowNames',true);
        
        if i==1&&j==1
            % only once: get reaction IDs from row names
            rxnIDs = tmpResTab.Properties.RowNames;
                        
            rxnNames = [model.rxnNames(contains(model.rxns,rxnIDs));
                {'uptake_reaction'}];
            rxnKEGGID = [model.rxnKEGGID(contains(model.rxns,rxnIDs));
                {''}];
            rxnKeggMaps = [model.rxnKeggMaps(contains(model.rxns,rxnIDs));
                {''}];
            subSystems = [model.subSystems(contains(model.rxns,rxnIDs)); {{'Medium'}}];
        end
        
        % convert both sampling tables to matrices
        fbaFluxes(:,c+1) = table2array(tmpResTab(:,1));
        fluxSamplingMat(:,colIdx+(n*c)) = table2array(tmpResTab(:,2:end));
        
        c=c+1;
        
    end
end; clear tmpResTab c colIdx

% set concentration values below 1E-13 to 0 (not treated as active in MILP)
concSamplingMat(concSamplingMat<1E-13) = 0;

% set flux values below 1E-15 to 0 (assumed numerical precision)
fluxSamplingMat(fluxSamplingMat<1E-15) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% two-way ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~~~~~~~~~~~~~~~~~~~~~ two-way anova for each gene ~~~~~~~~~~~~~~~~~~~~~ %

% Eg_i ~ carbonSource + uptakeFlux + carbonSource:uptakeFlux

% define factor names
factors = {'Source','Concentration','Interaction'};

% set significant level
alpha = 0.05;

% initialize matrix for P-values
pValues = ones(nGenes, 3);

% counter for all-zero concentration proteins
c=0;
for g=1:nGenes
    
    % create matrix for anova2 function
    y = reshape(concSamplingMat(g,:),...
        n*numel(concentrations),...
        numel(carbonSources));
    
    if sum(sum(y))==0
        c=c+1;
        fprintf('all zero concentrations for protein #%d\n', g)
    else
        pValues(g,:) = anova2(y, n, 'off');
    end
    
end; clear y g

pValues = pValues(~all(pValues==1,2),:);

fprintf('==> in total %d proteins had a concentration of zero across all conditions\n\n', c)

% adjust alpha for MHT correction (Bonferroni)
alpha = alpha / numel(pValues);

% ~~~ influence per subsystem ~~~ %

% remove orphan reactions
M = removeRxns(model,setdiff(model.rxns,rxnIDs));
M = removeRxns(M,findOrphanRxns(M));

% get list of subsystems
allSubSystems = vertcat(M.subSystems{:});

% unique subsystems in model
uniqueSubSystems = unique(allSubSystems);

% number of reactions per subsystem
nTotalSubsystems = cellfun(@(x)sum(ismember(allSubSystems, x)),...
    uniqueSubSystems);

% initialize matrix with influence of each factor on each subsystem
influencePerSubSystem = zeros(numel(uniqueSubSystems), size(pValues,2));

for i=1:size(pValues,2)
    
    % find significant p-values for the current coefficient
    idxSignificant = pValues(:,i)<alpha;
    
    
    if sum(idxSignificant)>0
        
        % find reactions associated to these genes
        tmpGenes = M.genes(idxSignificant);
        tmpRxns = findRxnsFromGenes(M, tmpGenes);
        tmpRxns = struct2cell(structfun(@(x)vertcat(x(:,1)), tmpRxns, 'un', 0));
        tmpRxns = unique(vertcat(tmpRxns{:}));
        
        % find subsystems associated to these reactions
        tmpSubSystems = [M.subSystems{contains(M.rxns, tmpRxns)}]';
        
        % count occurrences of each subsystem
        nPerSubsystem = cellfun(@(x)sum(ismember(tmpSubSystems, x)),...
            uniqueSubSystems);
        
        % scale to total numbers per subsystem
        influencePerSubSystem(:,i) = nPerSubsystem ./ nTotalSubsystems;
        
    end
end; clear tmpRxns tmpGenes dependentSubSystems ...
    pValues i idxSignificant idxNonSignificant M tmpIdx nPerSubSystem

% write ANOVA results to file
writetable(array2table(influencePerSubSystem, 'VariableNames', factors,...
    'RowNames', uniqueSubSystems),...
    [outDir, 'subsystems-anova-proteins.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');

% ~~~~~~~~~~~~~~~~~~~ two-way anova for each reaction ~~~~~~~~~~~~~~~~~~~ %
% remove blocked reactions from this analysis
blockedIdx = all(fluxSamplingMat<=0,2);
nzFluxMat = fluxSamplingMat(~blockedIdx,:);
nzSubSystems = subSystems(~blockedIdx);
nzRxnIDs = rxnIDs(~blockedIdx);

% set significance niveau
alpha = 0.05;

% initialize matrix of P-values
pValues = ones(size(nzFluxMat,1), 3); % v_i ~ carbonSource + uptakeFlux + carbonSource:uptakeFlux

for r=1:size(nzFluxMat,1)
    
   
    % create matrix for anova2 function
    y = reshape(nzFluxMat(r,:),...
        n*numel(concentrations),...
        numel(carbonSources));
     for i=1:size(y,2)
         for j=1:size(y,1)/n
             kstest(y((j-1)*n+1:j*n,i));
         end
     end
         
    pValues(r,:) = anova2(y, n, 'off');
    
end; clear y r blockedIdx

% adjust alpha for MHT correction (Bonferroni)
alpha = alpha / numel(pValues);

% ~~~ influence per subsystem ~~~ %

indRxns = {'ID','NAME','KEGGID','KEGGMAPS','TERM'};

% get list of subsystems
allSubSystems = vertcat(nzSubSystems{:});

% unique subsystems in model
uniqueSubSystems = unique(allSubSystems);

% number of reactions per subsystem
nTotalSubsystems = cellfun(@(x)sum(ismember(allSubSystems, x)),...
    uniqueSubSystems);

influencePerSubSystem = zeros(numel(uniqueSubSystems), size(pValues,2));

for i=1:size(pValues,2)
    % find significant p-values for the current coefficient
    idxSignificant = pValues(:,i)<alpha;
    
    if sum(idxSignificant)>0
        
        % find subsystems associated to these reactions
        tmpRxns = nzRxnIDs(idxSignificant);
        tmpSubSystems = [nzSubSystems{contains(nzRxnIDs, tmpRxns)}]';
        
        % count occurrences of each subsystem
        nPerSubsystem = cellfun(@(x)sum(ismember(tmpSubSystems, x)),...
            uniqueSubSystems);
        
        % scale to total numbers per subsystem
        influencePerSubSystem(:,i) = nPerSubsystem ./ nTotalSubsystems;
        
    end
end; clear tmpRxns tmpGenes nPerSubSystem nTotalSubsystems dependentSubSystems ...
    i alpha idxSignificant idxNonSignificant M tmpIdx nzSubSystems ...
    nzRxnIDs
% write anova results to file
writetable(array2table(influencePerSubSystem, 'VariableNames', factors,...
    'RowNames', uniqueSubSystems),...
    [outDir, 'subsystems-anova-reactions.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse concentration variability (CV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~ variance in concentration for every gene product over all conditions ~~~ %
cvConc = zeros(nGenes, 1);

for i=1:nGenes
    
    % total variation over all conditions
    tmpRow = concSamplingMat(i,:);
    
    cvConc(i) = std(tmpRow) / mean(tmpRow);
    
end; clear tmpRow

% ~~~ variation per subsystem ~~~ %
subSystemCV = zeros(nGenes,numel(uniqueSubSystems));

for i=1:numel(uniqueSubSystems)
    
    % find all genes associated with the current subsystem
    tmpRxns = findRxnsFromSubSystem(model, uniqueSubSystems{i});
    rxnIdx = findRxnIDs(model,tmpRxns);
    tmpGenes = findGenesFromRxns(model, tmpRxns);
    tmpGenes = unique(vertcat(tmpGenes{:}));
    
    if ~isempty(tmpGenes)
        % get all CV values for the genes found in the step before
        rowIdx=ismember(model.genes,tmpGenes);
        subSystemCV(rowIdx,i) = cvConc(rowIdx);
    end
end; clear tmpRxns tmpGenes rowIdx

colNames = regexprep(uniqueSubSystems, ' ', '_');
writetable(array2table(subSystemCV, 'VariableNames', colNames),...
    [outDir, 'protein-variances-per-subsystem.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', false, 'Delimiter', '\t');
clear colNames

% write whole set of abundance CV to file
writetable(array2table(cvConc,'VariableNames',{'CV'},'RowNames',model.genes),...
    [outDir,'abundance-CV.csv'], 'WriteRowNames', true)

% reactions associated to proteins with high variance
idxHighVarGenes = cvConc >= quantile(cvConc,0.9);
tmpRxns = findRxnsFromGenes(model,model.genes(idxHighVarGenes));
tmpRxns = struct2cell(structfun(@(x)vertcat(x(:,1)),tmpRxns,'un',0));
hvr = unique(vertcat(tmpRxns{:}));
clear tmpRxns idxHighVarGenes

tmpIdx=ismember(rxnIDs,hvr);
writetable(cell2table([rxnIDs(tmpIdx),rxnNames(tmpIdx),rxnKEGGID(tmpIdx),...
    rxnKeggMaps(tmpIdx),vertcat(subSystems{tmpIdx})],...
    'VariableNames', {'ID','NAME','KEGGID','KEGGMAP','SUBSYSTEM'}),...
    [outDir, 'high-variance-proteins-per-subsystem.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', false, 'Delimiter', '\t');
clear tmpIdx hvr

% reactions associated to proteins with low variance (remove proteins with
% all-zero concentrations)
idxLowVarGenes = cvConc < quantile(cvConc,0.1) & any(concSamplingMat,2);
tmpRxns = findRxnsFromGenes(model,model.genes(idxLowVarGenes));
tmpRxns = struct2cell(structfun(@(x)vertcat(x(:,1)),tmpRxns,'un',0));
lvr = unique(vertcat(tmpRxns{:}));
clear tmpRxns idxLowVarGenes

% do not consider blocked reactions
tmpIdx = ismember(rxnIDs,lvr) & any(fluxSamplingMat<=0,2);
writetable(cell2table([rxnIDs(tmpIdx),rxnNames(tmpIdx),rxnKEGGID(tmpIdx),...
    rxnKeggMaps(tmpIdx),vertcat(subSystems{tmpIdx})],...
    'VariableNames', {'ID','NAME','KEGGID','KEGGMAP','SUBSYSTEM'}),...
    [outDir, 'low-variance-proteins-per-subsystem.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', false, 'Delimiter', '\t');

clear tmpIdx lvr


% ~~~ variance in flux for every reaction over all conditions ~~~ %
cvFlux = zeros(nRxns, 1);

for i=1:nRxns
    
    % total variation over all conditions
    tmpRow = fluxSamplingMat(i,:);
    cvFlux(i) = std(tmpRow) / mean(tmpRow);
    
end; clear tmpRow
cvFlux(cvFlux<0) = 0;

% ~~~ variation per subsystem ~~~ %
subSystemCV = zeros(nRxns,numel(uniqueSubSystems));

for i=1:numel(uniqueSubSystems)
    
    % find reactions associated with current subsystem
    tmpRxns = findRxnsFromSubSystem(model, uniqueSubSystems{i});
    
    rowIdx = ismember(rxnIDs,tmpRxns);
    subSystemCV(rowIdx,i) = cvFlux(rowIdx);
    
end; clear tmpRxns rowIdx

colNames = regexprep(uniqueSubSystems, ' ', '_');
writetable(array2table(subSystemCV, 'VariableNames', colNames),...
    [outDir, 'reaction-variances-per-subsystem.csv'],...
    'WriteVariableNames', true, 'WriteRowNames', false, 'Delimiter', '\t');
clear colNames

% write whole set of flux CV to file
writetable(cell2table([model.rxns, model.rxnNames, vertcat(subSystems{:}),...
    model.rxnKEGGID, num2cell(cvFlux)],'VariableNames', {'ID','name','subsystem',...
    'KeggID','CV'}),...
    [outDir,'flux-CV.csv'])

% reactions with high cv
idxHighVarRxns = cvFlux >= quantile(cvFlux,0.9);
fName = [outDir, 'high-variance-reactions-per-subsystem.csv'];
if sum(idxHighVarRxns)>0
    writetable(cell2table([rxnIDs(idxHighVarRxns),rxnNames(idxHighVarRxns),...
        rxnKEGGID(idxHighVarRxns), rxnKeggMaps(idxHighVarRxns),...
        vertcat(subSystems{idxHighVarRxns})],...
        'VariableNames', {'ID','NAME','KEGGID','KEGGMAP','SUBSYSTEM'}),...
        fName,...
        'WriteVariableNames', true, 'WriteRowNames', false, 'Delimiter', '\t');
elseif isfile(fName)
    delete(fName)
end
clear idxHighVarRxns fName

% reactions with low variance (remove blocked reactions)
idxLowVarRxns = cvFlux < quantile(cvFlux,0.1) & any(fluxSamplingMat>0,2);
fName = [outDir, 'low-variance-reactions-per-subsystem.csv'];
if sum(idxLowVarRxns)>0
    writetable(cell2table([rxnIDs(idxLowVarRxns), rxnNames(idxLowVarRxns),...
        rxnKEGGID(idxLowVarRxns), rxnKeggMaps(idxLowVarRxns),...
        vertcat(subSystems{idxLowVarRxns})],...
        'VariableNames', {'ID','NAME','KEGGID','KEGGMAP','SUBSYSTEM'}),...
        fName,...
        'WriteVariableNames', true, 'WriteRowNames', false, 'Delimiter', '\t');
elseif isfile(fName)
    delete(fName)
end
clear idxLowVarRxns fName


% create a mapping between genes and reactions to be able to plot CV of
% protein concentrations against CV of reaction fluxes
rxnProtPairCounter = 0;
Z = zeros(nRxns,nGenes);
% non-zero entries indicate an association of protein and reaction, the value
% of the entry indicated the index of the reaction-protein pair
for i=1:nRxns
    geneIdx = cellfun(@str2double,regexp(model.rules{i},'\d+','match'));
    rxnProtPairCounter = rxnProtPairCounter+numel(geneIdx);
    Z(i,geneIdx) = rxnProtPairCounter-numel(geneIdx)+1:rxnProtPairCounter;
end

% dimensions #genes x #reactions
cvFluxConcMatch = nan(nGenes,nRxns);

% iterate over genes
for i=1:nGenes
    %     row: genes columns: corresponding reaction fluxes from cvFLux
    cvFluxConcMatch(i,Z(:,i)>0) = cvFlux(Z(:,i)>0);
end

writetable(array2table(cvFluxConcMatch),...
    [outDir,'rxn-prot-match-cv.csv'],...
    'WriteVariableNames',false,'WriteRowNames',false,'Delimiter','\t');

writetable(cell2table(vertcat(subSystems{:})),...
    [topDir,'/results/','subsystems.lst'],...
    'WriteVariableNames',false,'WriteRowNames',false,'FileType','text');

writetable(array2table(cvConc),[outDir,'prot-cv.csv'],...
    'WriteVariableNames',false,'WriteRowNames',false,'Delimiter','\t');
