function [bestSol] = IndiFinder(metScores, drugIndic)
%IndiFinder uses a genetic algorithm to identify the best predictors of a
%drug indication using metabolites scores for each sample
%
%   [bestSol] = IndiFinder(metScores, drugIndic)
%
%INPUTS
%   metScores   (m*x) matrix of the scores for each sample 
%               (m scores, x samples)
%
%   drugIndic   (1*x) vector of whether each sample has the given drug indication
%
%OUTPUTS
%   bestSol     vector (m*1) {0,1} containing the metabolites in the best
%               signature
%
% Originally coded by Jan Schellenberger as an implementation of the OptGene
% Algorithm as part of the COBRA Toolbox 2.0
%
% Modified to create IndiFinder by Daniel Zielinski, 07/2016

indicsLog = logical(drugIndic);
numIndics = sum(indicsLog);

global MaxPredictors

%Set the maximum number of predictors. The default behavior is up to 20
%samples containing the indication, use the number of positive samples as
%the number of maximum predictors. 
MaxPredictors = numIndics;
if MaxPredictors >20
    MaxPredictors = 20;
end

global HTABLE
HTABLE = java.util.Hashtable;

nmets = size(metScores,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS - set parameters here %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mutationRate = 1/nmets; % this should be changed as probably not a high enough rate
crossovermutationRate = mutationRate*.2;  % the rate of mutation after a crossover.  This value should probably be fairly low.  It is only there to ensure that not every member of the population ends up with the same genotype.
CrossoverFraction = .80;    % Percentage of offspring created by crossing over (as opposed to mutation). 0.7 - 0.8 were found to generate the highest mean, but this can be adjusted.
PopulationSize = [125 125 125 125]; % it was found that an increase beyond 125 individuals did not improve the results significantly.
Generations = 5000;    % maximum number of generations to perform
TimeLimit =  30*60;  % global time limit in seconds
StallTimeLimit = 3600*24*1;   % Stall time limit (terminate after this much time of not finding an improvement in fitness)
StallGenLimit =  Generations;       % terminate after this many generations of not finding an improvement
% PlotFcns =  {@gaplotscores, @gaplotbestf, @gaplotscorediversity, @gaplotstopping}; % what to plot.
PlotFcns =  {};
crossfun = @(a,b,c,d,e,f) crossoverCustom(a,b,c,d,e,f,crossovermutationRate);
MigrationFraction = .1;   % how many individuals migrate (.1 * 125 ~ 12 individuals).
MigrationInterval = 100;  % how often individuals migrate from one population to another.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InitialPopulation = [];


options = gaoptimset(                                       ...
    'PopulationType', 'bitstring',                          ...
    'CreationFcn', @lowmutationcreation,                    ...
    'MutationFcn', {@mutationUniformEqual, mutationRate},   ...
    'PopulationSize', PopulationSize,                       ...
    'StallTimeLimit', StallTimeLimit,                       ...
    'TimeLimit', TimeLimit,                                 ...
    'PlotFcns', PlotFcns,                                   ...
    'InitialPopulation', InitialPopulation,                 ...
    'CrossoverFraction', CrossoverFraction,                 ...
    'CrossoverFcn', crossfun,                               ...
    'StallGenLimit', StallGenLimit,                         ...
    'Generations', Generations,                             ...
    'TolFun', 1e-10,                                        ...
    'Vectorize', 'on',                                      ...
    'MigrationFraction', MigrationFraction,                 ...
    'MigrationInterval', MigrationInterval);

FitnessFunction = @(x) indicFitness(x,metScores,indicsLog);

gap.fitnessfcn  = FitnessFunction;
gap.nvars = nmets;
gap.options = options;

[x,FVAL,REASON,OUTPUT,population, scores] = ga(gap);

bestSol.points = x;
bestSol.val = FVAL;
bestSol.pop = population;
bestSol.allScores = scores;

return;


%% Creation Function
% generates initial warmup with much lower number of mutations (on average
% one mutation per
function [Population] = lowmutationcreation(GenomeLength,FitnessFcn,options)
totalPopulation = sum(options.PopulationSize);

global MaxPredictors;

%CREATE INITIAL POPULATION HERE

PopulationTemp = rand(totalPopulation,GenomeLength)*2-1;
PopulationTemp(abs(PopulationTemp)<(1-MaxPredictors/(2*GenomeLength))) = 0;
Population = sign(PopulationTemp);

return;

%% Mutation Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mutation function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mutationChildren = mutationUniformEqual(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate)

global MaxPredictors;

if(nargin < 8)
    mutationRate = 0.01; % default mutation rate
end
mutationChildren = zeros(length(parents),GenomeLength);
for i=1:length(parents)
    child = thisPopulation(parents(i),:);
    predictors = sum(abs(child));
    mutationPoints = find(rand(1,length(child)) < mutationRate);
    newMutations = rand(length(mutationPoints),1)*2-1;
    newMutations(abs(newMutations)<(1/3)) = 0;
    newMutations = sign(newMutations);
    child(mutationPoints) = newMutations;
    
    if MaxPredictors > 0
        while(sum(abs(child(:)))> MaxPredictors)
            ind2 = find(child);
            removeindex = ind2(randi(1,1,length(ind2))+1);
            child(removeindex) = 0;
        end
    end
    
    % with 50% chance, you will have fewer knockouts after mutation
    % than before.  This is to stop aquiring so many mutations.
    if rand > .5 && predictors > 1
        while(sum(abs(child(:)))>= predictors)
            ind2 = find(child);
            removeindex = ind2(randi(1,1,length(ind2))+1);
            child(removeindex) = 0;
        end
    end
    mutationChildren(i,:) = child;
end
return;

%% Crossover Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% crossover function %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xoverKids  = crossoverCustom(parents,options,GenomeLength,FitnessFcn,unused,thisPopulation, mutationRate)
nKids = length(parents)/2;
% Extract information about linear constraints
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

global MaxPredictors;

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;

% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    for j = 1:GenomeLength
        if(rand > 0.5)
            xoverKids(i,j) = thisPopulation(r1,j);
        else
            xoverKids(i,j) = thisPopulation(r2,j);
        end
    end
    if MaxPredictors>0
        while(sum(abs(xoverKids(i,:)))> MaxPredictors)
            ind2 = find(xoverKids(i,:));
            removeindex = ind2(randi(1,1,length(ind2))+1);
            xoverKids(i,removeindex) = 0;
        end
    end
end
% also apply mutations to crossover kids...
xoverKids = mutationUniformEqual(1:size(xoverKids,1) ,[],GenomeLength,[],[],[],xoverKids,mutationRate);
return;


function fitScore = indicFitness(pop,scores,indicsLogical)


popsize = size(pop,1);
fitScore = zeros(1,popsize);
T = length(find(indicsLogical));
F = length(find(~indicsLogical));


scoresTest = scores'*pop'; %each column is an individual, each row is a score for a sample for the side effect

%Score each individual as the inner product of their predictors and the
%sample scores
numPreds = sum(logical(pop),2);
allIndics = zeros(size(scoresTest));
for i = 1:size(pop,1)
    curSides = sortrows([scoresTest(:,i),indicsLogical'],1);
    allIndics(:,i) = curSides(:,2);
end

%Calculate the AUC
indexMat = tril(ones(size(allIndics,1)));
TPs = +allIndics'*indexMat;
FPs = +(~allIndics)'*indexMat;
FNs = T-TPs;
TNs = F-FPs;
TPRs = TPs ./ (TPs+FNs);
FPRs = FPs ./ (FPs+TNs);

mutPen = log10(numPreds+1); %add one to prevent infinity error
%Assign a mutation penalty to remove neutral predictors
for i = 1:length(fitScore)
    fitScore(i) = trapz(FPRs(i,:),TPRs(i,:))*.9999^mutPen(i);
end

