function [bestmod_num_1, idxPredictors_1, termLabels_1, allData1, allData_table_d1, allData_table2_1, ...
    labels, idxOutcome, nparamfolds, inmod_min_allx, cvfits] ...
    = driver_function(ADHD, idxPredictors, idxOutcome, apply2test,k,thresh,lin)

allData=ADHD.allData;

allData1=allData;
allData_table_d1=ADHD.allData_table;
allData_table2_1=ADHD.allData_table2;
train_data=ADHD.train_data;
test_data=ADHD.test_data;

if lin==1
zsc = true; % standardize
end
if lin==0,
zsc = false; % standardize
end

do1se = false;
randomizeOutcome=false;

data = ADHD.train_data;
ID=ADHD.train_ID;
tdata = ADHD.test_data;
tID = ADHD.test_ID;

n=size(data,1);
nt=size(tdata,1);
nparamfolds=20;
nv=floor(n/k);
cvfits = cell(k,1);
Rsqs = zeros(k,2);
modstrs = cell(k,2);
lammin=.01;
nlam=100;

if(zsc),data=zscore(data);tdata=zscore(tdata);end
termLabels = ADHD.labels(idxPredictors);
labels=ADHD.labels;

%% Regress and tune parameters within nested k-fold

if lin==1,

%diary('output.txt');
for i=1:k
    
    traindata = data([1:nv*(i-1) nv*i+1:end],:);
    trainID = ID([1:nv*(i-1) nv*i+1:end]);
    valdata = data(nv*(i-1)+1:nv*i,:);
    valID = ID(nv*(i-1)+1:nv*i);
    
    % Set up design matrix with all terms
    design = traindata(:,idxPredictors);
    vdesign = valdata(:,idxPredictors);
    
    % Fit with elastic net
    [ cvfit,modstr,modstr1se,inmod_min,inmod_1se ] = fitADHD_eNet_FINAL( design,traindata(:,idxOutcome),termLabels,nparamfolds,nlam,lammin,idxOutcome,lin);
    cvfits{i}=cvfit;
   
    %NEW - save for all vars if in model
    inmod_min_all{i}=inmod_min;
    inmod_min_allx(i,:)=inmod_min;
    
    % track # models in which each variable appears
    if(i==1), modcount=double(inmod_min); else modcount = modcount + double(inmod_min); end
    if(i==1), modcount1se=double(inmod_1se); else modcount1se = modcount1se + double(inmod_1se); end
  
    
    % print model
    fprintf(1,'\n\n****************************** FOLD %d *****************************************\n\n',i);
    fprintf(1,'****************\n');
    fprintf(1,'\n\n* Model with min MSE *\n\n');
    fprintf(1,'****************\n\n');
    disp(modstr);
    
    % Predict on validation set
    [ AICc,Rsq ] = fitADHD_test( cvfit,vdesign,valdata(:,idxOutcome),thresh,'lambda_min');
    Rsqs(i,1) = Rsq;
    modstrs{i,1} = modstr;
    
    if(do1se)
        fprintf(1,'\n\n****************\n');
        fprintf(1,'\n\n* Model within 1 SE of min MSE *\n\n');
        fprintf(1,'****************\n\n');
        disp(modstr1se);
        
        % Predict on validation set
        [ AICc,Rsq ] = fitADHD_test( cvfit,vdesign,valdata(:,idxOutcome),thresh,'lambda_1se' );
        Rsqs(i,2) = Rsq;
        modstrs{i,2} = modstr1se;
    end
    
end
    
disp(['Mean validation R^2 across folds: ' num2str(mean(Rsqs(:,1))) ' +/- ' num2str(std(Rsqs(:,1)))]);
end

if lin==0,
    %diary('output.txt');
for i=1:k
    
    traindata = data([1:nv*(i-1) nv*i+1:end],:);
    trainID = ID([1:nv*(i-1) nv*i+1:end]);
    valdata = data(nv*(i-1)+1:nv*i,:);
    valID = ID(nv*(i-1)+1:nv*i);
    
    % Set up design matrix with all terms
    design = traindata(:,idxPredictors);
    vdesign = valdata(:,idxPredictors);
    
    % Fit with elastic net
    [ cvfit,modstr,modstr1se,inmod_min,inmod_1se ] = fitADHD_eNet_FINAL( design,traindata(:,idxOutcome),termLabels,nparamfolds,nlam,lammin,idxOutcome,lin);
    cvfits{i}=cvfit;
    
    %NEW - save for all vars if in model
    inmod_min_all{i}=inmod_min;
    inmod_min_allx(i,:)=inmod_min;
    
    % track # models in which each variable appears
    if(i==1), modcount=double(inmod_min); else modcount = modcount + double(inmod_min); end
    if(i==1), modcount1se=double(inmod_1se); else modcount1se = modcount1se + double(inmod_1se); end
  
% print model
    fprintf(1,'\n\n****************************** FOLD %d *****************************************\n\n',i);
    fprintf(1,'****************\n');
    fprintf(1,'\n\n* Model with min deviance *\n\n');
    fprintf(1,'****************\n\n');
    disp(modstr);
    
    % Predict on validation set
    [ AICc,Rsq ] = fitADHD_test( cvfit,vdesign,valdata(:,idxOutcome),thresh,'lambda_min');
    Rsqs(i,1) = Rsq;
    modstrs{i,1} = modstr;
    
    if(do1se)
        fprintf(1,'\n\n****************\n');
        fprintf(1,'\n\n* Model within 1 SE of min deviance *\n\n');
        fprintf(1,'****************\n\n');
        disp(modstr1se);
        
        % Predict on validation set
        [ AICc,Rsq ] = fitADHD_test( cvfit,vdesign,valdata(:,idxOutcome),thresh,'lambda_1se' );
        Rsqs(i,2) = Rsq;
        modstrs{i,2} = modstr1se;
    end
    
  end

disp(['Mean validation AUC across folds: ' num2str(mean(Rsqs(:,1))) ' +/- ' num2str(std(Rsqs(:,1)))]);
  end
%% Make parts to build a new model using covs from best model

%pick model for which rsquare is highest
[~,bestmod_num_1]= max(Rsqs(:,1));

%make an idxBrain/idxPredictor variable for only the variables in inmod_min
%for highest rsquare model
t=(inmod_min_all(1,bestmod_num_1));
idxPredictors_1=idxPredictors(t{1,1});

termLabels_1 = ADHD.labels(idxPredictors_1);

%%apply to test
if(apply2test)
    
    design = data(:,idxPredictors);
    tdesign = tdata(:,idxPredictors);
    
    % print model
    fprintf(1,'\n\n****************************** TEST *****************************************\n\n');
    fprintf(1,'******** final model: best from k-fold ********\n');
    fprintf(1,'****************\n\n');
    
    
    fprintf(1,'* Model with minimum MSE *\n\n');
    [maxRsq,idxmod] = max(Rsqs(:,1));
    cvfit = cvfits{idxmod};
    modstr = modstrs{idxmod,1};
    
    disp(modstr);
    
    % Predict on test set
    fitADHD_test( cvfit,tdesign,tdata(:,idxOutcome),thresh,'lambda_min' );
    
    if(do1se)
        fprintf(1,'\n\n* Model within 1 SE of min MSE *\n\n');
        [maxRsq,idxmod] = max(Rsqs(:,2));
        cvfit = cvfits{idxmod};
        modstr1se = modstrs{idxmod,2};
        
        disp(modstr1se);
        
        % Predict on test set
        fitADHD_test( cvfit,tdesign,tdata(:,idxOutcome),thresh,'lambda_1se' );
    end
    
%     %% Predict on test set, using selected vars, training on whole train set
    
    % restrict to selected vars
    svars=logical(modcount);%>k-4;
    design = design(:,svars); tdesign = tdesign(:,svars); termLabels = termLabels(svars);
    
    % Fit with elastic net
    [ cvfit,modstr,modstr1se,inmod_min,inmod_1se ] = fitADHD_eNet_FINAL( design,data(:,idxOutcome),termLabels,nparamfolds,nlam,lammin,idxOutcome,lin);
    % print model
    fprintf(1,'\n\n****************************** TEST *****************************************\n\n');
    fprintf(1,'******** final model: trained on whole training set ********\n');
    fprintf(1,'****************\n\n');
    
    
    fprintf(1,'* Model with minimum MSE *\n\n');
    disp(modstr);
    
    % Predict on test set
    fitADHD_test( cvfit,tdesign,tdata(:,idxOutcome) );
    
    if(do1se)
        fprintf(1,'\n\n* Model within 1 SE of min MSE *\n\n');
        disp(modstr1se);
        
        % Predict on test set
        fitADHD_test( cvfit,tdesign,tdata(:,idxOutcome),'lambda_1se' );
    end
    
end

%diary OFF