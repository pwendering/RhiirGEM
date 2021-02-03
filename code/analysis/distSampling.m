function samples = distSampling(model, n, minLimit, maxLimit)
%% samples = distSampling(model, n, minLimit, maxLimit)
% Sample n flux distributions by minimizing the distance of v to a random
% flux vector v* (contains values between minimal and maximal values for
% each reaction determined by FVA)
% Input:
%   struct model:               metabolic model structure
%   double n:                   desired number of sampled flux
%                               distributions
%   double minLimit:            minimal flux values for each reaction
%   double maxLimit:            maximal flux values for each reaction
% Output:
%   double samples:             matrix containg the sampled flux
%                               distributions (#reactions x #samples)

rng('default')

% initialize variables
nRxns = size(model.S,2);
nMets = size(model.S,1);
idxObj = model.c==1;
samples = zeros(nRxns,n);

if nargin<3
    minLimit=0;
    maxLimit=max(model.ub);
end

if numel(minLimit)==1
    minLimit = minLimit*ones(nRxns,1);
end

if numel(maxLimit)==1
    maxLimit = maxLimit*ones(nRxns,1);
end

%% set up the linear optimization program
problem.Aeq = model.S;
problem.beq = model.b;
problem.Aineq = [];
problem.bineq = [];
problem.lb = model.lb;
problem.ub = model.ub;

% solve FBA and set biomass to optimal value
s = optimizeCbModel(model);
maxGrowthRate = s.x(idxObj);
problem.lb(idxObj) = maxGrowthRate;
problem.ub(idxObj) = maxGrowthRate;
clear s

% add constraints for sampling with objective:
% min |v* - v| transformed to min sum( theta+ + theta-)
% -v - theta+ + theta- = -v*
thetaMatrix = [
    -eye(nRxns),...    v
    -eye(nRxns),...    theta+
    eye(nRxns)...      theta-
];

% inequality constraints LHS
problem.Aineq = thetaMatrix;

% adjust size of equality matrix
problem.Aeq = [problem.Aeq, zeros(nMets,2*nRxns)];

% define objective
problem.f = [zeros(nRxns,1); ones(2*nRxns,1)];

% append lower bounds for thetas
problem.lb = [problem.lb; zeros(2*nRxns,1)]; 

% append upper bounds for thetas
problem.ub = [problem.ub; ones(2*nRxns,1)]; 


%% sample n times
parfor i=1:n
    
    % random sample of fluxes between given limits, e.g. from FVA
    randFluxes = (maxLimit-minLimit) .* rand(nRxns,1) + minLimit;
    randFluxes(idxObj) = maxGrowthRate;
    
    % change the linear objective
    tmpProblem = problem;
    tmpProblem.bineq = -randFluxes;
    
    % solve the LP
    s = cplexlp(tmpProblem);
    
    if ~isempty(s)
        samples(:,i) = s(1:nRxns);
    end
end

end