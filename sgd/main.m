%% Tissue-level SGD
%%% HAM
[tD,mO,gE,lR] = getSGDinRefModels('Recon3','HAM','woFluxes');
[tD,mO,gE,lR] = getSGDinRefModels('HMR3765','HAM','woFluxes');
[tD,mO,gE,lR] = getSGDinRefModels('GBM','HAM','woFluxes');

%%% FBS
[tD,mO,gE,lR] = getSGDinRefModels('Recon3','FBS','woFluxes');
[tD,mO,gE,lR] = getSGDinRefModels('HMR3765','FBS','woFluxes');
[tD,mO,gE,lR] = getSGDinRefModels('GBM','FBS','woFluxes');