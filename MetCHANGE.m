function [objDM] = MetCHANGE(model,exprPVals,solver)
%MetCHANGE solves a series of linear programming problems to calculate the
%consistency of optimal metabolite production for each metabolite with gene
%expression data. These scores then must be compared between control and
%treated data using a standard score to obtain MetCHANGE perturbation
%scores (not included). Must have the COBRA Toolbox and a supported
%linear programming solver package in order to run this program.
%
%   [objDM] = MetCHANGE(model,exprPVals,solver)
%
%INPUTS
%   model   (the following fields are required - others can be supplied)
%       S       Stoichiometric matrix
%       b       Right hand side = dx/dt
%       c       Objective coefficients
%       lb      Lower bounds
%       ub      Upper bounds
%       mets    List of metabolites
%       rxns    List of reactions
%
%   exprPVals   presence/absence p-values for each reaction (not gene)
%
%   solver      Name of linear programming solver package loaded
%         (From COBRA Toolbox):
%         Currently allowed LP solvers:
%           lindo_new       Lindo API >v2.0
%           lindo_old       Lindo API <v2.0
%           glpk            GLPK solver with Matlab mex interface (glpkmex)
%           lp_solve        lp_solve with Matlab API
%           tomlab_cplex    CPLEX accessed through Tomlab environment (default)
%           cplex_direct    CPLEX accessed direct to Tomlab cplex.m. This gives
%                   the user more control of solver parameters. e.g.
%                   minimising the Euclidean norm of the internal flux to
%                   get rid of net flux around loops
%           mosek           Mosek LP solver with Matlab API (using linprog.m included in Mosek
%                   package)
%           gurobi          Gurobi accessed through Matlab mex interface (Gurobi mex)
%           matlab          Matlab's own linprog.m (currently unsupported, may not
%                   work on COBRA-type LP problems)
%           mps             Outputs a MPS matrix string. Does not solve LP problem
%
%OUTPUTS
%   objDM       matrix (m*x) containing objective scores for each sample
%               (m scores, x samples)
%
% Written by Daniel Zielinski and Monica Mo
% Date 6/24/2012


%Set solver tolerance for off reactions
tol = 10^-7;

%Get reaction IDs and set up objective
expression = exprPVals;
nConditions=size(expression,2);

%Make model irreversible
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);

%Option to set certain uptakes to 0
removeUptakes = {};

%Make pvalue irreversible
nIrrevRxns = size(irrev2rev,1);
expressionIrrev = zeros(nIrrevRxns,nConditions);
for i=1:nIrrevRxns
    expressionIrrev(i,:) = expression(irrev2rev(i,1),:);
end

%LP 1 from Supplementary Figure 1
%Get optimal productions for each metabolite
maxObjectiveDM = zeros(size(modelIrrev.mets,1),1);
for i = 1:length(modelIrrev.mets)
    changeCobraSolver(solver,'LP');
    modelIrrevTemp = modelIrrev;
    if ~isempty(removeUptakes)
        modelIrrevTemp = changeRxnBounds(modelIrrevTemp,removeUptakes,0,'u');
    end
    [modelIrrevTemp,rxnIDexists] = addReaction(modelIrrevTemp,'DMSpec_cur',modelIrrevTemp.mets(i),-1,true,0,1000,0,'','');
    if isempty(rxnIDexists)
        modelIrrevTemp = changeObjective(modelIrrevTemp,modelIrrevTemp.rxns(end),1);
    else
        modelIrrevTemp = changeObjective(modelIrrevTemp,modelIrrevTemp.rxns(rxnIDexists),1);
    end
    %find max objective
    FBAsolution = optimizeCbModel(modelIrrevTemp);
    if (abs(FBAsolution.f)<tol)
        %Need upper bounds on EX_*_b and DM_ reactions to be 1
        exchRxns = [intersect(strmatch('EX_',modelIrrevTemp.rxns),find(~cellfun(@isempty,strfind(modelIrrevTemp.rxns,'_b'))));strmatch('DM_',modelIrrevTemp.rxns)];
        if~isempty(rxnIDexists)&&~isempty(find(exchRxns==rxnIDexists, 1))
            exchRxns = setxor(rxnIDexists,exchRxns);
        end
        modelIrrevTemp = changeRxnBounds(modelIrrevTemp,modelIrrevTemp.rxns(exchRxns),1,'u');
        modelIrrevTemp = changeRxnBounds(modelIrrevTemp,modelIrrevTemp.rxns(exchRxns),0,'l');
        FBAsolution = optimizeCbModel(modelIrrevTemp);
    end
    if (FBAsolution.stat ~= 1)
        maxObjectiveDM(i)=0;
    end
    maxObjectiveDM(i)=FBAsolution.f;
end

%Determine which metabolites cannot be produced, to avoid wasting time
%calculating their solutions
toRemove = [find(maxObjectiveDM==1000);find(abs(maxObjectiveDM)<tol)];
toRemove = sort(toRemove);


%LP 2 from Supplementary Figure 1
%Get consistencies of optimal production with gene expression presence/absence p values
objDM = zeros(length(modelIrrev.mets),nConditions);
for i=1:length(modelIrrev.mets)
    if isempty(find(toRemove==i, 1))
        curObj = maxObjectiveDM(i);
        modelIrrevTemp = modelIrrev;
        if ~isempty(removeUptakes)
            modelIrrevTemp = changeRxnBounds(modelIrrevTemp,removeUptakes,0,'u');
        end
        [modelIrrevTemp,rxnIDexists] = addReaction(modelIrrevTemp,'DMSpec_cur',modelIrrevTemp.mets(i),-1,true,0,1000,0,'','');
        if isempty(rxnIDexists)
            modelIrrevTemp = changeObjective(modelIrrevTemp,modelIrrevTemp.rxns(end),1);
        else
            modelIrrevTemp = changeObjective(modelIrrevTemp,modelIrrevTemp.rxns(rxnIDexists),1);
        end
        sol = [];
        for j=1:nConditions
            changeCobraSolver(solver,'LP');
            model2gimme=modelIrrevTemp;
            if isempty(rxnIDexists)
                model2gimme.c=[expressionIrrev(:,j);0];
                model2gimme.lb(end) = curObj;
                model2gimme.ub(end) = curObj;
            else
                model2gimme.c= expressionIrrev(:,j);
                model2gimme.lb(rxnIDexists) = curObj;
                model2gimme.ub(rxnIDexists) = curObj;
            end
            gimmeSolution = optimizeCbModel(model2gimme,'min');
            if (abs(gimmeSolution.f)<tol)
                exchRxns = [intersect(strmatch('EX_',model2gimme.rxns),find(~cellfun(@isempty,strfind(model2gimme.rxns,'_b'))));strmatch('DM_',model2gimme.rxns)];
                if ~isempty(rxnIDexists)&&~isempty(find(exchRxns==rxnIDexists, 1))
                    exchRxns = setxor(rxnIDexists,exchRxns);
                end
                model2gimme = changeRxnBounds(model2gimme,model2gimme.rxns(exchRxns),1,'u');
                model2gimme = changeRxnBounds(model2gimme,model2gimme.rxns(exchRxns),0,'l');
                gimmeSolution = optimizeCbModel(model2gimme,'min');
            end
            if (gimmeSolution.stat ==1 )
                sol(j)=gimmeSolution.f;
            else
                sol(j) = 0;
            end
        end
        objDM(i,:) = sol;
    end
end