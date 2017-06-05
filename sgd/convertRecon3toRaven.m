clear all

%Load src and model
addpath(genpath('src/'))
load 'models/Recon3d.mat'
model = modelRecon3model;

%Remove redundant reacton
I     = find(strcmp(model.rxns,'EX_M01874[e]'));
model = removeRxns(model,model.rxns(I(2:end)));

%Make it Raven compatible
model.metNames = model.mets;

%%% Add comps
model = addMetcomp(model);
tempmets  = model.mets;
tempNames = model.metNames;
for j = 1:length(model.mets)
    curmet = strsplit(model.mets{j}, '['); %remove trailing [c]
    model.mets{j} = curmet{1};
end
model.metNames = model.mets;

%Close model - this means to take all exchange rxns and provide a fake
%compartment in which mets are added to block exchange of mets in exchange
%rxns. E.g. 13_cis_oretn[n] => becomes 13_cis_oretn[n] => 13_cis_oretn[b]
%where [b] is the fake compartment. 13_cis_oretn[b] is flagged as
%unconstrained in the model structure.
model.comps             = [model.comps;'b'];
[I,J]                   = getExchangeRxns(model);
metsToAdd.mets          = strcat(I,'_fake');
metsToAdd.metNames      = metsToAdd.mets;
metsToAdd.compartments  = 'b';
metsToAdd.unconstrained = true(numel(I),1);
newmodel                = addMets(model,metsToAdd);
newmodel.S(:,J)         = [model.S(:,J);eye(numel(J))];

%Check that all exchange reactions are stoichiometrically closed
printModel(newmodel,getExchangeRxns(newmodel))
printModel(simplifyModel(newmodel),getExchangeRxns(newmodel))
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';
selExc = (find( full((sum(abs(newmodel.S)==1,1) ==1) & (sum(newmodel.S~=0) == 1))))';
model = newmodel;

%%Add id
model.id = 'Recon3D';
model.description = 'Recon3D in RAVEN format';

%Check tasks
taskReport = checkTasks(model, 'data/recon3D_Tasks.xlsx', true);

%Save
save('models/Recon3d_toRaven','model')