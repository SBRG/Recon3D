
function [HMR2Reconmets_mapped,HMR2Reconmets_unmapped] = mapMetsHMR2Recon(metList)
    load models/HMRmodel.mat
    HMRmodel = HMRmodel;
    load models/2016_07_13_Recon3.mat
    recon3d  = modelRecon3;
    metMap = zeros(length(HMRmodel.mets), 1);
    for i = 1:length(HMRmodel.mets)
        if ~isnan(HMRmodel.chebi{i})
           currentCompartment = HMRmodel.comps{HMRmodel.metComps(i)};
           if strcmp(currentCompartment,'s')
               currentCompartment = 'e';
           end
           currentChebi = HMRmodel.chebi{i};

           ReconMet = find(ismember(recon3d.metCHEBIID, currentChebi));
           for j = 1:length(ReconMet)
               ReconCompartment = strsplit(recon3d.mets{ReconMet(j)}, '[');
               ReconCompartment = strsplit(ReconCompartment{2},']');
               ReconCompartment = ReconCompartment{1};
               if strcmp(ReconCompartment, currentCompartment)
                    metMap(i) = ReconMet(j);
                    break
                end
           end          
        end
    end
    sum(sign(metMap))

    metNames = modifyMetNames(HMRmodel);
    
    %Map names
    metListindinHMR = [];
    for i = 1:length(metList)
        metListindinHMR = [metListindinHMR; find(ismember(metNames,metList(i)))];
    end
    HMR2Reconmets          = metMap(metListindinHMR);
    HMR2Reconmets_mapped   = recon3d.mets(HMR2Reconmets(HMR2Reconmets~=0));
    HMRmetlist_mapped      = metList(HMR2Reconmets~=0);
    HMR2Reconmets_unmapped = metList(HMR2Reconmets==0);
    
    %print mapping
    fprintf('%d%s\n%s\t%s\n',length(HMR2Reconmets_mapped),' metabolites were mapped from HMR to Recon:','HMR','Recon')
    for i = 1:length(HMR2Reconmets_mapped)
           fprintf('%s\t%s\n', HMRmetlist_mapped{i}, HMR2Reconmets_mapped{i})
    end
    
    fprintf('%d%s\n',length(HMR2Reconmets_unmapped),' metabolites were not mapped from HMR to Recon:')
    for i = 1:length(HMR2Reconmets_unmapped)
           fprintf('%s\n', HMR2Reconmets_unmapped{i})
    end
end
