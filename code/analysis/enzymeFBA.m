function [solution,problem] = enzymeFBA(model, kcats, mw, C, objIdx, verbose)
%% solution = enzymeFBA(model, kcats, mw, C, objIdx(opt), verbose(opt))
% MILP that imposes upper bounds on reaction fluxes by v_max =
% kcat * [E]
% The concentration of every protein is integrated with the use of GPR
% rules taking that minimum for AND rules and the sum for OR rules.
% The total mass of enzymes is constrained by a given value C.
% Input:
%       struct model:           metabolic model, osenseStr is used to
%                               decide on minimization or maximization
%       double kcats:           maximum turnover for every reaction
%       double mw:              molecular weights for every protein
%       double C:               upper bound on protein mass
%       double objIdx:          (optional) index of the objective to 
%                               optimize
%       logicel verbose:        (optional) whether or not to print output
% Output:
%       struct solution:        output from cplexlp
%               x           solution found by the optimization function.
%                           if exitflag > 0, then x is a solution;
%                           otherwise, x is the value of the optimization 
%                           routine when it terminated prematurely.
%
%               fval        optimal growth rate [x(objIdx)]
%
%               exitflag 	integer specifying the reason the optimization 
%                           algorithm terminated
%
%               output      structure containing information about the 
%                           optimization.
%                           The fields of the structure are:
%
%                   iterations          Number of iterations
%                   algorithm           Optimization algorithm used
%                   message             Exit message
%                   time                Execution time of the algorithm
%                   cplexstatus         Status code of the solution
%                   cplexstatusstring   Status string of the solution
% 
%       struct problem:        MILP which was run to obtain the solution x
%                              Fields:
% 
%                   f           objective
%                   Aeq         equality constraints
%                   beq         target values for equality constraints
%                   Aineq       inequality constraints
%                   bineq       target values for inequality constraints
%                   lb          lower bounds
%                   ub          upper bounds
%                   ctype       variable types (continuous/binary)

if nargin < 5 || isempty(objIdx)
    objIdx = find(model.c);
end

if nargin < 6
    verbose = false;
end

if ~isfield(model, 'osenseStr')
    warning('model does not have the field ''osenseStr'', adding max objective sense')
    model.osenseStr = 'max';
end

nGenes = numel(model.genes);
nRxns = size(model.S,2);
nMets = size(model.S,1);

% make mw a row vector
mw = reshape(mw,1,numel(mw));

%% add constraints for reaction upper bounds and enzyme concentrations
[Eineq, bineq, Eeq, beq] = getEnzymeConstraints(model, kcats);

% number of reaction-protein pairs
tmp=cellfun(@(x)regexp(x,'\d+','match'),model.rules,'un',0);
tmp=[tmp{:}];
nRxnProtPairs = numel(tmp);
clear tmp

% number of additional alternative complexes
nANDComps = size(Eineq,2)-2*nGenes-nRxns-nRxnProtPairs;

% define upper bounds; reaction flux cannot be larger than the largest kcat
% as the enzyme concentration will always be below 1 mmol gDW-1
upperBounds = [max(kcats)*ones(nRxns,1);ones(2*nGenes+nRxnProtPairs+nANDComps,1)];

% find boundary reacions (upper bounds should remain the same)
excIdx = find(findExcRxns(model));
upperBounds(excIdx) = model.ub(excIdx);

% find NGAM reaction and keep upper bound 
ngamIdx = contains(model.rxnNames,'ngam');
if ~any(ngamIdx)
    error('no NGAM found using query ''ngam''')
end
upperBounds(ngamIdx)=model.ub(ngamIdx);

% objective
f = zeros(size(Eineq,2),1);
if isequal(model.osenseStr,'max')
    f(objIdx) = -1;
elseif isequal(model.osenseStr,'min')
    f(objIdx) = 1;
else
    error('No objective sense found')
end

%% set up the MILP

% equality constraints
% LHS
problem.Aeq = sparse([
    model.S,... v
    zeros(nMets,nGenes),... E_g
    zeros(nMets,nRxnProtPairs),... reaction-protein pairs
    zeros(nMets,nGenes),... y
    zeros(nMets,nANDComps);... AND-components
    Eeq... enzyme equality constraints (i.e. promiscuity constraints)
]);

% RHS
problem.beq = [
    zeros(nMets,1);... steady state
    beq... promiscuity constraints
];

% inequality constraints

% LHS
problem.Aineq = sparse([
    zeros(1,nRxns),... v
    mw,... E_g
    zeros(1,nRxnProtPairs),... reaction-protein pairs
    zeros(1,nGenes),... y
    zeros(1,nANDComps);... AND-components
    Eineq... enzyme inequality constraints
]);

% RHS
problem.bineq = [
    C;... protein mass fraction as upper bound for enzyme mass fractions
    bineq... enzyme inequality constraints
];

problem.lb = [model.lb; zeros(2*nGenes+nRxnProtPairs+nANDComps,1)];
problem.ub = upperBounds;
problem.f = f;

% define variable types (y is binary variable)
problem.ctype = [
    repmat('C',1,nRxns+nGenes+nRxnProtPairs),... v,E_G,protRxnPairs
    repmat('B',1,nGenes),... y
    repmat('C',1,nANDComps)... AND-components
];

% pass options to the CPLEX solver
problem.options = cplexoptimset('cplex');
% set integer tolerance to 1E-10 (default: 1E-5)
problem.options.mip.tolerances.integrality = 1E-10;
% set feasibility and optimality tolerances to 1E-10 (default: 1E-6)
problem.options.simplex.tolerances.feasibility = 1E-10;
problem.options.simplex.tolerances.optimality =  1E-10;

[solution.x, ~, solution.exitflag, solution.output] = cplexmilp(problem);

if verbose
    disp(solution.output.message)
    if ~isempty(solution.x)
        % index ranges for variables
        rxnIdx=1:nRxns;
        proteinIdx=nRxns+1:nRxns+nGenes;
        prPairIdx=nRxns+nGenes+1:nRxns+nGenes+nRxnProtPairs;
        yIdx=nRxns+nGenes+nRxnProtPairs+1:size(Eineq,2)-nANDComps;
        
        fprintf('number of active reactions: %d\n', sum(solution.x(rxnIdx)>0)) % v > 0
        fprintf('number of translated proteins: sum(y): %d\n',...
            sum(solution.x(yIdx))) % sum y
        fprintf('sum E-R pairs - sum E_g: %d\n',...
            sum(solution.x(prPairIdx))-sum(solution.x(proteinIdx))) % precision of promiscuity constraint
        fprintf('total protein mass: %.4f (C=%.4f)\n',...
            sum(solution.x(proteinIdx).*mw'),... % check if mass-constraint is fulfilled
            C)
        fprintf('objective value: %.5f\n', solution.x(objIdx)) % growth value
    end
end

end

