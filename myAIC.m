function AIC = myAIC(RSS,n,k,corrected)
%	function AIC = myAIC(RSS,n,k,corrected)
%   Returns Akaike Information Criterion, corrected if corrected=true (default).

if(nargin<4), corrected=true; end

AIC = n*log(RSS/n) + 2*k;
if(corrected), AIC = AIC + 2*k*(k+1)/(n-k-1); end


end

