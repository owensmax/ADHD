function [ cvfit,modstr_min,modstr_1se,inmod_min,inmod_1se ] = fitADHD_eNet_RR( design,comp,termLabels,nparamfolds,nlam,lammin,idxOutcome,lin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Defaults

if(nargin<6), lammin=.01; end
if(nargin<5), nlam=100; end
if(nargin<4), nparamfolds=10; end
if(nargin<7), idxOutcome=4;end %added for logistic regression

%% Fit with elastic net

opt = glmnetSet;
opt.alpha=0.01;
opt.nlambda=nlam;%1000;
opt.lambda_min=lammin;%.1;
alphas = 1:-0.05:0.05;
%%for linear regression
if lin==1, cvfit =cvALglmnet(design,comp,'gaussian',opt,'default',nparamfolds,alphas);end
%for logistic regression
if lin==0, cvfit = cvALglmnet(design,comp,'multinomial',opt,'default',nparamfolds,alphas);end
% cvglmnetPlot(cvfit)
% tlocx = mean(get(gca,'xlim'));
% tlocy = max(get(gca,'ylim')) - .1*range(get(gca,'ylim'));
% text(tlocx,tlocy,['\alpha = ' num2str(cvfit.alpha)]);
% drawnow
% set(gcf,'color','white')
%end

%% Make model strings

% for min MSE
coeffs=cvglmnetCoef(cvfit,'lambda_min');
inmod_min = (coeffs(2:end)~=0)';
if lin==1,%linear regression
modstr_min = sprintf('cbcl_adhd_r = %3.2g',coeffs(1));end
if lin==0%logistic regression
modstr_min = sprintf('log-odd(outcome) = %3.2g',coeffs(1));end
for j = 2:length(coeffs)
    if(coeffs(j))
        modstr_min = [modstr_min sprintf(' + %3.2g*%s',coeffs(j),termLabels{j-1})];
    end
end

% for within 1se of min MSEs
coeffs=cvglmnetCoef(cvfit,'lambda_1se');
inmod_1se = (coeffs(2:end)~=0)';
if lin==1,%linear regression
modstr_1se = sprintf('cbcl_adhd_r = %3.2g',coeffs(1));end
if lin==0,%logistic regression
modstr_1se = sprintf('log-odd(outcome) = %3.2g',coeffs(1));end
for j = 2:length(coeffs)
    if(coeffs(j))
        modstr_1se = [modstr_1se sprintf(' + %3.2g*%s',coeffs(j),termLabels{j-1})];
    end
end

end

