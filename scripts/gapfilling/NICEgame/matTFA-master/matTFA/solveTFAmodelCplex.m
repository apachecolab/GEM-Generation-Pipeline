function sol = solveTFAmodelCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)
% Solve a model using specific solver settings. More details in
% changeToCPLEX_WithOptions
%
% INPUTS
%   tModel: a TFA-ready model (has a .A matrix)
%   TimeInSec: timelimit for the solver. 
%   manualScalingFactor: manual scaling factor for the solver
%   mipTolInt: Integer tolerance of the solver
%   emphPar: Solver emphasis (trade-offs between speed, feasibility, 
%       optimality, and moving bounds in MIP - see https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/MIPEmphasis.html)
%   feasTol: solver tolerance for feasibility (error on constraints)
%   mipDisplay: verbosity of the MIP info display
%   CPXPARAMdisp: Turn on/off the cplex problem setup display
% 
%% Changelog
% 2017/04/26 - Modified by Pierre on Georgios Fengos's base, to incorporate
% Vikash's Gurobi hooks in a more global fashion, in a similar way COBRA
% solvers are handled.
% Calls the global parameter TFA_MILP_SOLVER, and check if it set to
% something. If not, default to CPLEX using Fengos's code.
% In the long run, we should rename thins function to remove CPLEX from its
% name.
%
% 2026/04/14 - Modified by Sudharshan Ravi to enable execution on Apple
% Silicon (MACA64). CPLEX is x86-only and cannot load on ARM64 MATLAB.
% Default solver switched from cplex_direct to gurobi_direct. Full Gurobi
% parameter block added to x_solveGurobi to mirror CPLEX LCSB defaults:
% IntFeasTol=1e-9, FeasibilityTol=1e-9, NumericFocus=3 (matches CPLEX
% emphPar=1 extreme caution), ScaleFlag=0 (no scaling, matches scalPar=-1),
% TimeLimit forwarded when set. Warning message updated to include Gurobi
% status string. Missing semicolons fixed in x_solveGurobi.
%
global TFA_MILP_SOLVER

if ~exist('TFA_MILP_SOLVER','var') || isempty(TFA_MILP_SOLVER)
    TFA_MILP_SOLVER = 'gurobi_direct'; % default changed from cplex_direct: CPLEX is x86-only, cannot load on MACA64
end
solver = TFA_MILP_SOLVER;
% solver = 'gurobi_direct';
if ~exist('manualScalingFactor','var') || isempty(manualScalingFactor)
    manualScalingFactor = [];
end
if ~exist('mipTolInt','var') || isempty(mipTolInt)
    mipTolInt = [];
end
if ~exist('emphPar','var') || isempty(emphPar)
    emphPar = [];
end
if ~exist('feasTol','var') || isempty(feasTol)
    feasTol = [];
end
if ~exist('scalPar','var') || isempty(scalPar)
    scalPar = [];
end
if ~exist('TimeInSec','var') || isempty(TimeInSec)
    TimeInSec = [];
end
if ~exist('mipDisplay','var') || isempty(mipDisplay)
    mipDisplay = [];
end
if ~exist('CPXPARAMdisp','var') || isempty(CPXPARAMdisp)
    CPXPARAMdisp = [];
end

switch solver
    %% Case CPLEX
    case 'cplex_direct'
        if isempty(which('cplex.p'))
            error('You need to add CPLEX to the Matlab-path!!')
        end
        sol = x_solveCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp);
    %% Case GUROBI
    case 'gurobi_direct'
        if isempty(which('gurobi'))
            error('You need to add Gurobi to the Matlab-path!!')
        end
        sol = x_solveGurobi(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay);
end
end

%% Private function for CPLEX solve
function sol = x_solveCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)

% this function solves a TFBA problem using CPLEX

% if cplex is installed, and in the path
if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

% Convert problem to cplex
cplex = changeToCPLEX_WithOptions(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp);

% Optimize the problem
try
    CplexSol = cplex.solve();
    if isfield(cplex.Solution,'x')
        x = cplex.Solution.x;
        if ~isempty(x)
            sol.x = cplex.Solution.x;
            sol.val = cplex.Solution.objval;
            sol.cplexSolStatus = cplex.Solution.status;
        else
            sol.x = [];
            sol.val = [];
            disp('Empty solution');
            warning('Cplex returned an empty solution!')
            sol.cplexSolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.val = [];
        disp('No field cplex.Solution.x');
        warning('The solver does not return a solution!')
        sol.cplexSolStatus = 'No field cplex.Solution.x';
    end
catch
    sol.x = NaN;
    sol.val = NaN;
    sol.cplexSolStatus = 'Solver crashed';
end

delete(cplex)
end

%% Private function for GUROBI solve
function sol = x_solveGurobi(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay)
    
    num_constr = length(tModel.constraintType);
    num_vars = length(tModel.vartypes);

    contypes = '';
    vtypes = '';

    % convert contypes and vtypes into the right format
    for i=1:num_constr
        contypes = strcat(contypes,tModel.constraintType{i,1});
    end
    
    for i=1:num_vars
        vtypes = strcat(vtypes,tModel.vartypes{i,1});
    end
    
    gmodel.A=tModel.A;
    gmodel.obj=tModel.f;
    gmodel.lb=tModel.var_lb;
    gmodel.ub=tModel.var_ub;
    gmodel.rhs=tModel.rhs;
    gmodel.sense=contypes;
    gmodel.vtype=vtypes;
  
    gmodel.varnames=tModel.varNames;
    
    if tModel.objtype==-1
      gmodel.modelsense='max';
    elseif tModel.objtype==1
        gmodel.modelsense='min';
    else
        error(['No objective type specified ' ...
            '(model.objtype should be in {-1,1})']);
    end

    params.OutputFlag = 0;
    params.IntFeasTol     = 1e-9;
    params.FeasibilityTol = 1e-9;
    params.NumericFocus   = 3;
    params.ScaleFlag      = 0;
    if ~isempty(TimeInSec)
        params.TimeLimit = TimeInSec;
    end
    try
        result=gurobi(gmodel,params);
        if isfield(result,'x')
            x = result.x;
            x(find(abs(x) < 1E-9))=0;
        else
            warning('The solver does not return a solution! Gurobi status: %s', result.status);
            result.x=[];
            result.objval=[];
        end
    catch
        result.status='0';
        x=NaN;
        result.x=NaN;
        result.objval=NaN;
    end
    
    sol.x=result.x;
    sol.val=result.objval;
    sol.status=result.status;

end