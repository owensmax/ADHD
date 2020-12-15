clear all;

files=["/home/max/Documents/ABCD_ADHD/smri_win_RR_standardcovs.csv", "standcov_ADHDdata_mri";...
    "/home/max/Documents/ABCD_ADHD/sst_win_RR_standardcovs.csv", "standcov_ADHDdata_sst";...
    "/home/max/Documents/ABCD_ADHD/nb_win_RR_standardcovs.csv", "standcov_ADHDdata_nb";...
    "/home/max/Documents/ABCD_ADHD/mid_win_RR_standardcovs.csv", "standcov_ADHDdata_mid"]

for i=1:size(files,1)

    allData=readtable(files(i,1));%Read in Data From Rstudio

    [envir]=preanalysis_function_FINAL(allData,i);
    
    load envir; %load envir saved in function
    save(files(i,2)); %name final of environment

end
