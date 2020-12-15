function [ AICc,Rsq ] = fitADHD_test( cvfit,tdesign,tcomp,thresh,mod )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if(nargin<4), mod = 'lambda_min'; end

nt=length(tcomp);

if(strcmp(cvfit.name,'Binomial Deviance'))
    probs = cvglmnetPredict(cvfit,tdesign,mod,'response');
    [FPR, TPR, T, AUC, OPTROCPT] = perfcurve(logical(tcomp),probs,true);
    nP = sum(tcomp); nN = length(tcomp)-nP;
    %disp([nP nN])
    p = auc2p(AUC,nP,nN,thresh);
    optFPR = OPTROCPT(1); optTPR = OPTROCPT(2);
    TP = optTPR*nP; FP = optFPR*nN;
    PPV = TP / (TP+FP);
    F1 = 2*(PPV*optTPR)/(PPV+optTPR);
    coeffs=cvglmnetCoef(cvfit,mod);
    DF = nnz(coeffs)+1;
    fprintf('\n')
    disp(['Stats for model from elastic net applied to test data, ' num2str(DF) ' DOF:'])
    disp([ 'AUC = ' num2str(AUC)])
    %disp([ 'Rbar^2 = ' num2str(1 - SSresid*(nt-1)/(nt-DF-1)/SStotal)])
    disp(['p = ' num2str(p)]);    
    AICc = p; Rsq = AUC;
else
    SStotal = sum((tcomp - mean(tcomp)).^2);
    pred = cvglmnetPredict(cvfit,tdesign,mod);
    SSresid = sum((tcomp - pred).^2);
    coeffs=cvglmnetCoef(cvfit,mod);
    DF = nnz(coeffs)+1;
    AICc = myAIC(SSresid,nt,DF);
    Rsq = 1 - SSresid/SStotal;
    fprintf('\n')
    disp(['Stats for model from elastic net applied to test data, ' num2str(DF) ' DOF:'])
    disp([ 'R^2 = ' num2str(Rsq)])
    %disp([ 'Rbar^2 = ' num2str(1 - SSresid*(nt-1)/(nt-DF-1)/SStotal)])
    disp(['AICc = ' num2str(AICc)]);

end
end

