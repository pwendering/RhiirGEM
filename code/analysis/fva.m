function [minFlux,maxFlux] = fva(model)
%% [minFLux, maxFlux] = fva(model)
% Flux Variability Analysis (FVA)
% Input:
%       struct model:           metabolic model
% Output:
%       double minFLux:         vector of minimal fluxes
%       double maxFlux:         vector of maximal fluxes

nRxns = size(model.S,2);

% initialize output variables
minFlux = zeros(nRxns,1);
maxFlux = zeros(nRxns,1);

problem.Aeq=model.S;
problem.Aineq=[];
problem.beq=model.b;
problem.bineq=[];
problem.lb=model.lb;
problem.ub=model.ub;
clear model

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