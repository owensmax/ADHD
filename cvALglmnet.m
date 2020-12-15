function cvfit = cvALglmnet(x,y,family,options,type,nfolds,alphas,foldid,plots)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%% Set default values
if nargin < 3 || isempty(family)
    family = 'gaussian';
end
if nargin < 4
    options = [];    
end
if nargin < 5 || isempty(type)
    type = 'default';
end
if nargin < 6 || isempty(nfolds)
    nfolds = 10;
end
if nargin < 7
    alphas = 1:-0.01:0.01;
end
nalpha = length(alphas);
if nargin < 8
    foldid = [];
end
if nargin < 9
    plots = false;
end

alphamat = repmat(alphas,options.nlambda,1);
lambdamat = zeros(options.nlambda,nalpha);
cvm = zeros(options.nlambda,nalpha);
cvsd = zeros(options.nlambda,nalpha);
cvfits = cell(nalpha,1);

options.alpha = alphas(1);
cvfitbest = cvglmnet(x,y,family,options,type,nfolds,foldid,false,true);
cvfitbest.alpha = alphas(1);
cvfits{1} = cvfitbest;
if(length(cvfitbest.cvm)==options.nlambda)
    cvm(:,1) = cvfitbest.cvm;
    cvsd(:,1) = cvfitbest.cvsd;
    lambdamat(:,1) = cvfitbest.lambda;
else
    wrongn = length(cvfitbest.cvm);
    cvm(1:wrongn,1) = cvfitbest.cvm; cvm(wrongn+1:end,1)=NaN;
    cvsd(1:wrongn,1) = cvfitbest.cvsd; cvsd(wrongn+1:end,1)=NaN;
    lambdamat(1:wrongn,1) = cvfitbest.lambda; lambdamat(wrongn+1:end,1)=NaN;
end
foldid = cvfitbest.foldid;

for i=2:nalpha
    
    options.alpha = alphas(i);
    cvfitnew = cvglmnet(x,y,family,options,type,nfolds,foldid,false,true);
    cvfitnew.alpha = alphas(i);
    cvfits{i} = cvfitnew;
    if(length(cvfitnew.cvm)==options.nlambda)
        cvm(:,i) = cvfitnew.cvm;
        cvsd(:,i) = cvfitnew.cvsd;
        lambdamat(:,i) = cvfitnew.lambda;
    else
        wrongn = length(cvfitnew.cvm);
        cvm(1:wrongn,i) = cvfitnew.cvm; cvm(wrongn+1:end,i)=NaN;
        cvsd(1:wrongn,i) = cvfitnew.cvsd; cvsd(wrongn+1:end,i)=NaN;
        lambdamat(1:wrongn,i) = cvfitnew.lambda; lambdamat(wrongn+1:end,i)=NaN;
    end
    if( min(cvfitnew.cvm) < min(cvfitbest.cvm) ), cvfitbest = cvfitnew; end
    %nnz(cvglmnetCoef(cvfitnew,'lambda_min')) <= nnz(cvglmnetCoef(cvfitbest,'lambda_min')) && ...
    
end

if(cvfitbest.nzero(cvfitbest.lambda==cvfitbest.lambda_1se))

idx = find(cvfitbest.lambda==cvfitbest.lambda_min);
target = cvfitbest.cvm(idx)+ cvfitbest.cvsd(idx);
idx = find(cvfitbest.cvm < target,1);

if(~isempty(idx))

targetDF = cvfitbest.nzero(idx);

for i=1:nalpha
   
    cvfitcand = cvfits{i};
    idx = find(cvfitcand.cvm < target,1);
    if(~isempty(idx) && cvfitcand.nzero(idx))
        if(cvfitcand.nzero(idx) < targetDF)
            cvfitbest = cvfitcand;
            cvfitbest.lambda_1se = cvfitbest.lambda(idx);
            targetDF = cvfitbest.nzero(idx);
        end
    end
    
end

end

end

cvfit = cvfitbest;

if(plots)
    
plot3_errorbars_surf(alphamat(:), lambdamat(:), cvm(:), cvsd(:));   
    
end

end

