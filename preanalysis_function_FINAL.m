function [envir] = driver_function(allData_p,i)

%clearvars
close all;
truei=i;

%load in template
load ADHDdata;
allData=allData_p;


%make SID
allID=allData(:,1);
allID_t=allID;
allData.id_redcap=nominal(allData.id_redcap);

allData_table=[allData(:,2:10) allData(:,11:width(allData))];

%Make Datasets into Numeric Data
for i = 2:width(allData_table)
    var = allData_table.(i);
    if iscellstr(var)
       allData_table.(i) = str2double(var);
    end
end

allData=table2array(allData_table);

allData_table2=[allID allData_table];

% Check for missing data

TF = isnan(allData);
TFsums=sum(TF);
missing_data=sum(sum(TF(:,:)))

%make labels
labels=allData_table.Properties.VariableNames;
labels=labels(:,1:width(allData_table));

%make ids
idxBrain=(2:width(allData_table));
%idxKSADS_adhd=1;
idxCBCL_att=1;

%Make Training and Test Datasets

train_IDt=cell2table(train_ID);
test_IDt=cell2table(test_ID);

test_IDt.Properties.VariableNames(1) = "id_redcap";
test_IDt=innerjoin(test_IDt,allID_t);
test_ID=table2cell(test_IDt);

train_IDt=setdiff(allID_t,test_IDt);
train_ID=table2cell(train_IDt);

test_data=innerjoin(test_IDt,allData_table2);
test_data=test_data(:,2:width(allData_table2));
test_data=table2array(test_data(:,1:width(allData_table)));

train_data=innerjoin(train_IDt,allData_table2);
train_data=train_data(:,2:width(allData_table2));
train_data=table2array(train_data(:,1:width(allData_table)));

%remove extraneous variables
clear idxCBCL missing_data allData_table3 testidx trainidx test_IDt train_IDt TF TFsums fs_data ans demos ndvars test_ID_mid var rdVals seVals siVals allData_table_array;
i=truei;
save envir;
envir=load("envir.mat");

