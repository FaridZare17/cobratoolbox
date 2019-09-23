% The COBRAToolbox: testRemoveZeroFluxRxns.m
%
% Purpose:
%     - testremoveZeroFluxRxns.m tests the functionality of
%     removeZeroFluxRxns and checks to see whether removing reactions
%     by this function changes in phenotype of the microorganism
%
% Authors:
%     -  Farid Zare 9/15/2019
%

% save the current path
currentDir = pwd;

% initialize the test
fileDir = fileparts(which('testRemoveZeroFluxRxns'));
cd(fileDir);

% load model
model = getDistributedModel('ecoli_core_model.mat');

% define the solver packages to be used to run this test
solverPkgs = { 'mosek', 'gurobi', 'tomlab_cplex', 'glpk'};

for k = 1:length(solverPkgs)
    
    % change the COBRA solver (LP)
    solverOK = changeCobraSolver(solverPkgs{k}, 'LP', 0);
    
    if solverOK == 1
        fprintf('   Testing "removeZeroFluxRxns" function, using %s ... ', solverPkgs{k});
        
        % set the tolerance
        tol = getCobraSolverParams('LP','feasTol');
        
        % Remove reactions using function
        modelOut = removeZeroFluxRxns(model);
        
        % testing if the objectvie value is unchanged
        sol = optimizeCbModel(model);
        sol1 = optimizeCbModel(modelOut);
        assert(abs(sol.f-sol1.f)<tol);
        
        % testing if all the reactions with value > tol still exist in the
        % model
        id = abs(sol.x) > tol ; 
        assert(sum(ismember(modelOut.rxns, model.rxns(id))) == length(model.rxns(id))) ;
        
        % testing the phenotype, no significant change in any fluxes in the
        % curated model
        rxnsID = findRxnIDs(model, modelOut.rxns);
        solFlux = sol.x(rxnsID) ;
        assert(all(abs(solFlux - sol1.x) < tol));
        
        % check if any reaction is removed to see if deleting reactions works properly
        assert(length(modelOut.rxns) ~= length(model.rxns));
        
        % check if the functions omits the Fructose exchange reaction which
        % contains no flux
        assert(~any(contains(modelOut.rxns,{'EX_fru(e)'})));
        
%         % check if the function works correctly for the optional inputs
        
%         % reactions MALt2_2 and EX_gln_L(e) cannot have any flux and must be
%         % removed, and other 2 reactions must remain in the model
%         rxnNames = {'MALt2_2', 'EX_gln_L(e)', 'GAPD', 'ACALD'};
%         modelOut1 = removeZeroFluxRxns(model, rxnNames);
%         assert(~any(ismember(modelOut1.rxns,{'MALt2_2', 'EX_gln_L(e)'})));
%         assert(any(ismember(modelOut1.rxns,{'GAPD'})));
%         assert(any(ismember(modelOut1.rxns,{'ACALD'})));
%         
%         % testing if all other reactions still exist in the model
%         assert(length(modelOut1.rxns) + 2 == length(model.rxns));
        
        % check if the upper and lower bounds are set correctly
        % the maximum uptake for oxygen is 60 and minimun is 0
        o2ID = findRxnIDs(modelOut, 'EX_o2(e)');
        assert(modelOut.lb(o2ID) + 60 < tol);
        assert(modelOut.ub(o2ID) < tol);
        
        % output a success massage
        fprintf('Done.\n');
    end
end

% change the directory
cd(currentDir)
