function [minFlux,maxFlux] = fva(model,f)
%% [minFLux, maxFlux] = fva(model)
% Flux Variability Analysis (FVA)
% Input:
%       struct model:           metabolic model
%       double f:               percentage of optimal growth at which FVA
%                               will be performed (default: 1)
% Output:
%       double minFlux:         vector of minimal fluxes
%       double maxFlux:         vector of maximal fluxes

if nargin < 2
    f = 1;
elseif f > 1
    f = f / 100;
end

nRxns = size(model.S,2);
bioIdx = find(model.c==1);

% initialize output variables
minFlux = zeros(nRxns,1);
maxFlux = zeros(nRxns,1);

problem.Aeq=model.S;
problem.Aineq=[];
problem.beq=model.b;
problem.bineq=[];
problem.lb=model.lb;
problem.ub=model.ub;
problem.f=-model.c;
clear model

% pass options to the CPLEX solver
problem.options = cplexoptimset('cplex');
% set integer tolerance to 1E-10 (default: 1E-5)
problem.options.mip.tolerances.integrality = 1E-10;
% set feasibility and optimality tolerances to 1E-10 (default: 1E-6)
problem.options.simplex.tolerances.feasibility = 1E-10;
problem.options.simplex.tolerances.optimality =  1E-10;

x = cplexlp(problem);
problem.lb(bioIdx) = f*x(bioIdx);

parfor i=1:nRxns
    
    tmpProblem=problem;
    
    f = zeros(nRxns,1);
    
    % minimization
    f(i) = 1;
    tmpProblem.f = f;
    s = cplexlp(tmpProblem);
    minFlux(i) = s(i);
    
    % maximization
    f(i) = -1;
    tmpProblem.f = f;
    s = cplexlp(tmpProblem);
    maxFlux(i) = s(i);
    
end

end