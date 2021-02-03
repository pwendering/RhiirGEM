function [eConcMat, rxnFluxMat] = enzymeSampling(model, minLimit, maxLimit, kcats, mw, C, n)
%% enzymeSampling(model, minLimit, maxLimit, kcats, mw, C, n)
% minimize the distance between a vector of random concentrations and the
% predicted abundances
% Input:
%   struct model:           metabolic model
%   double minLimit:        lower abundance limit for each individual 
%                           protein or common lower limit for all proteins 
%                           (will be used to find random abundances within 
%                           the lower and upper limit)
%   double maxLimit:        similar to lowerLimit
%   double kcats:           turnover numbers for every reaction [h^-1]
%   double mw:              molecular weights for each protein [g mmol^-1]
%   double C:               total protein content [g/gDW]
%   double n:               desired number of sampled distributions
% Output:
%   double eConcMat:        matrix containing protein abundances for each
%                           sampled point
%   double rxnFluxMat:      matrix of reaction fluxes for each sampled
%                           point

% initialize variables
nRxns = size(model.S,2);
nGenes = numel(model.genes);
objIdx = model.c==1;
eConcMat = zeros(nGenes,n);
rxnFluxMat = zeros(nRxns,n);

if numel(minLimit)==1
    minLimit = minLimit*ones(nGenes,1);
end

if numel(maxLimit)==1
    maxLimit = maxLimit*ones(nGenes,1);
end

%% find optimal growth rate and reference concentrations
[s,problem] = enzymeFBA(model, kcats, mw, C, [], false);
problem.lb(objIdx) = 0.99*s.x(objIdx);
clear s

% number of reaction-protein pairs
tmp = cellfun(@(x)regexp(x,'\d+','match'),model.rules,'un',0);
tmp = [tmp{:}];
nRPP = numel(tmp);
clear tmp

% number of AND-components
nANDComps = size(problem.Aineq,2)-nRxns-2*nGenes-nRPP;

% E_g + theta+ - theta- = E_g*
thetaMatrix = [
    zeros(nGenes,nRxns),...     v
    eye(nGenes),...             E_g
    zeros(nGenes,nRPP),...      RPP
    zeros(nGenes,nGenes),...    y
    zeros(nGenes,nANDComps),... AND-components
    eye(nGenes),...             theta+
    -eye(nGenes)...             theta-
];

% add matrix for thetas to the equality matrix
problem.Aeq = sparse([
    [problem.Aeq, zeros(size(problem.Aeq,1),2*nGenes)];...
    thetaMatrix
]);

% add an empty matrix to the right of the inequality matrix
problem.Aineq = sparse([
    problem.Aineq, zeros(size(problem.Aineq,1),2*nGenes)...
]);

% append lower and upper bounds for delta+ and delta-
problem.lb = [
    problem.lb;...
    zeros(2*nGenes,1)...
];

problem.ub = [
    problem.ub;...
    ones(2*nGenes,1)...
];

% define the objective
problem.f = [
    zeros(size(problem.f));...  v, E_g, RPP, y, AND-Components
    ones(nGenes,1);...          theta+
    ones(nGenes,1)...           theta-
];

problem.ctype = [problem.ctype, repmat('C',1,2*nGenes)];

%% sample n times
parfor i=1:n
    
    tmpProblem = problem;
    
    % random sample of concentrations
    randConc = (maxLimit-minLimit) .* rand(nGenes,1) + minLimit;
    
    % RHS vector for equality constraints:
    tmpProblem.beq = [tmpProblem.beq; randConc];
    x = cplexmilp(tmpProblem);

    if ~isempty(x)
        % store result in the output matrices
        rxnFluxMat(:,i) = x(1:nRxns);
        eConcMat(:,i) = x(nRxns+1:nRxns+nGenes);
    end
end

end