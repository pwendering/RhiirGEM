function [minConc, maxConc, minFlux, maxFlux] = enzymeFVA(model, kcats, mw, C, includeFlux)
%% [minC, maxC] = enzymeFVA(model, kcats, mw, C)
% Variability analysis for the enzyme concentrations using enzymeFBA
% The concentration range of every protein is determined while ensuring
% optimal biomass production
% Input:
%       struct model:           metabolic model
%       double kcats:           maximum turnover for every reaction
%       double mw:              molecular weights for every protein
%       double C:               upper bound on protein mass
%       logical includeFlux:    whether to perform FVA additionally
% Output:
%       double minConc:         vector of minimal protein concentrations
%       double maxConc:         vector of maximal protein concentrations

nGenes = numel(model.genes);
nRxns = size(model.S,2);
objIdx = find(model.c==1);

if nargin < 5
    includeFlux = 0;
    minFlux = nan(nRxns,1);
    maxFlux = minFlux;
end

% get optimal growth value
[s,problem] = enzymeFBA(model, kcats, mw, C, [], false);
model.lb(objIdx) = 0.99*s.x(objIdx);
clear s

% initialize output variables
maxConc = zeros(nGenes,1);
minConc = zeros(nGenes,1);

% reset objective
problem.f(:)=0;

parfor i=(nRxns+1):(nRxns+nGenes)
    tmpProblem = problem;
    
    % maximization
    tmpProblem.f(i) = -1;
    s = cplexmilp(tmpProblem);
    maxConc(i-nRxns) = s(i);
    
    % minimization
    tmpProblem.f(i) = 1;
    s = cplexmilp(tmpProblem);
    minConc(i-nRxns) = s(i);
end

if includeFlux
    parfor i=1:nRxns
        tmpProblem = problem;
        
        % maximization
        tmpProblem.f(i) = -1;
        s = cplexmilp(tmpProblem);
        maxFlux(i) = s(i);
        
        % minimization
        tmpProblem.f(i) = 1;
        s = cplexmilp(tmpProblem);
        minFlux(i) = s(i);
    end
end
end

