%pick which runs

clear *
load ADHDdata_mri

smri=false;
nback=false;
sst=true;
mid=true;

apply2test = true;

lin=true;
log=false;

thresh=.5;

k=5;
%% Smri
if (smri),
    
ADHD_mri = load("standcov_ADHDdata_mri.mat");

%set predictors to be used
idxPredictors = ADHD_mri.idxBrain;

if (lin), idxOutcome = ADHD_mri.idxCBCL_att;end
if (log), idxOutcome = ADHD_mri.idxKSADS_adhd;end

%initiate the driver function, modified slightly from the original driver script from Nick
[bestmod_num_1, idxPredictors_1, termLabels_1, allData1, allData_table_1, allData_table2_1, ...
    labels_1, idxOutcome, nparamfolds,inmod_min_all_1,cvfits_1] = driver_function_FINAL(ADHD_mri,idxPredictors,idxOutcome,...
    apply2test,k,thresh(1),lin);

end;
%% Nback
if (nback),
ADHD_nb = load("standcov_ADHDdata_nb.mat");

idxPredictors = ADHD_nb.idxBrain;

if (lin), idxOutcome = ADHD_nb.idxCBCL_att;end
if (log), idxOutcome = ADHD_nb.idxKSADS_adhd;end

[bestmod_num_2, idxPredictors_2, termLabels_2, allData4, allData_table_2, allData_table2_2, ...
    labels_2, idxOutcome, nparamfolds, inmod_min_all_2] = driver_function_FINAL(ADHD_nb,idxPredictors,idxOutcome,...
    apply2test,k,thresh,lin);

end
%% SST
if (sst),
ADHD_sst = load("standcov_ADHDdata_sst.mat");

idxPredictors = ADHD_sst.idxBrain;

if (lin), idxOutcome = ADHD_sst.idxCBCL_att;end
if (log), idxOutcome = ADHD_sst.idxKSADS_adhd;end

[bestmod_num_3, idxPredictors_3, termLabels_3, allData5, allData_table_3, allData_table2_3, ...
    labels_3, idxOutcome, nparamfolds, inmod_min_all_3] = driver_function_FINAL(ADHD_sst,idxPredictors,idxOutcome,...
    apply2test,k,thresh,lin);

end

%% MID
if (mid),
ADHD_mid = load("standcov_ADHDdata_mid.mat");

idxPredictors = ADHD_mid.idxBrain;

if (lin), idxOutcome = ADHD_mid.idxCBCL_att;end
if (log), idxOutcome = ADHD_mid.idxKSADS_adhd;end

[bestmod_num_4, idxPredictors_4, termLabels_4, allData6, allData_table_4, allData_table2_4,...
    labels_4, idxOutcome, nparamfolds, inmod_min_all_4] = driver_function_FINAL(ADHD_mid,idxPredictors,idxOutcome,...
    apply2test,k,thresh,lin);

end
