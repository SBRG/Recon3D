function mapMetsRecon2HMR(metList)    
    load ../models/HMRmodel.mat
    HMRmodel = HMRmodel;
    load ../models/2016_07_13_Recon3.mat
    recon3d  = modelRecon3;
    metMap = zeros(length(recon3d.mets), 1);
    for i = 1:length(recon3d.mets)
        if ~isnan(recon3d.metCHEBIID{i})
           currentCompartment = strsplit(recon3d.mets{i}, '[');
           currentCompartment = strsplit(currentCompartment{2},']');
           currentCompartment = currentCompartment{1};
           currentChebi = recon3d.metCHEBIID{i};

           HMRMet = find(ismember(HMRmodel.chebi, currentChebi));
           for j = 1:length(HMRMet)
                HMRCompartment = HMRmodel.metComps(HMRMet(j));
                HMRCompartment = HMRmodel.comps{HMRCompartment};
                if strcmp(HMRCompartment,'s')
                    HMRCompartment = 'e';
                end

                if strcmp(HMRCompartment, currentCompartment)
                    metMap(i) = HMRMet(j);
                    break
                end
           end          
        end
    end
    sum(sign(metMap))

    metNames = modifyMetNames(HMRmodel);
    
    %Map names
    metListindinRecon = [];
    for i = 1:length(metList)
        metListindinRecon = [metListindinRecon; find(ismember(recon3d.mets,metList(i)))];
    end
    HMR2Reconmets          = metMap(metListindinRecon);
    HMR2Reconmets_mapped   = metNames(HMR2Reconmets(HMR2Reconmets~=0));
    Reconmetlist_mapped    = metList(HMR2Reconmets~=0);
    HMR2Reconmets_unmapped = metList(HMR2Reconmets==0);
    
    %print mapping
    fprintf('%d%s\n%s\t%s\n',length(HMR2Reconmets_mapped),' metabolites were mapped from HMR to Recon:','Recon','HMR')
    for i = 1:length(HMR2Reconmets_mapped)
           fprintf('%s\t%s\n', Reconmetlist_mapped{i}, HMR2Reconmets_mapped{i})
    end
    
    fprintf('%d%s\n',length(HMR2Reconmets_unmapped),' metabolites were not mapped from HMR to Recon:')
    for i = 1:length(HMR2Reconmets_unmapped)
           fprintf('%s\n', HMR2Reconmets_unmapped{i})
    end
end

function metNames = modifyMetNames(modelMod)
% modifyMetNames
% Adds the compartment to the metabolite names e.g. 'glucose' in
% compartment s becomes 'glucose[s]'.
%
%   modelMod         a model structure
%   metNames         metabolite names with compartment indications
%
%   Avlant Nilsson, 2016-05-16
%
    metNames = modelMod.metNames;
    for i=1:length(metNames)
        metNames{i} = [metNames{i} '[' modelMod.compNames{modelMod.metComps(i)} ']'];
    end
end

