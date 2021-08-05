function gurobiProblem = cplex2gurobi(problem)
% converts a cplex-formatted optimization problem to gurobi format
% customized for this paper, probably not applicable for the general case
problem.obj = abs(problem.f);
problem.A = [problem.Aeq; problem.Aineq];
problem.sense = [repmat('=',1,size(problem.Aeq,1)) repmat('<',1,size(problem.Aineq,1))];
problem.rhs = [problem.beq; problem.bineq];
problem.vtype = problem.ctype;
problem.modelsense = 'max';
gurobiProblem = rmfield(problem,{'Aeq','Aineq','beq','bineq','f','ctype','options'});
end