function K = mysgwt_func(kid)
%MYSGWT_FUNC return wavelet function Ks.
%   
%   K includes the kernel function g and normalization factor.
%
%   K := (\sum_k \lambda_k) * ( g(s lambda_i)/s ) /  (sum_j g(s
%   lambda_j)/s)
%      = (\sum_k \lambda_k) * ( g(s lambda_i)) /  (sum_j g(s
%   lambda_j))
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2015/11/01 19:45:32 (CST) $
    switch kid
        case 1
            K = @kernel1;
        case 2
            K = @kernel2;
    end
            
end


function [kvals, dkvals] = kernel1(s, lambdaa)
%   Kernel funciton with normalization.
%   K = @(s, x) x.* exp(-sx);
%   ldd := K(s lambda)/s
    ldd = lambdaa.* exp(-s.*lambdaa);
    sumld = sum(lambdaa); % L1 norm of lambdaa
    sumldd = sum(ldd);
    kvals = sumld.* ldd/sumldd;
    
% dksx = dK/ds [K(s, x)]    
    dksx = @(s, ld) -ld.^2.* exp(-s*ld); 
    T1 = dksx(s,lambdaa);
    dkvals = sumld* ( T1./sumldd + ldd*sumabs(T1)/sumldd^2);
end

function [kvals, dkvals] = kernel2(s, lambdaa)
%   Kernel funciton without normalization.
%   K = @(s, x) x.* exp(-sx);
%   ldd := K(s lambda)/s
    ldd = lambdaa.* exp(-s.*lambdaa);
    kvals = ldd;
    
% dksx = dK/ds [K(s, x)]    
    dksx = @(s, ld) -ld.^2.* exp(-s*ld); 
    dkvals = dksx(s,lambdaa);
end

