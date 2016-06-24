function [A, s, ghistory, fhistory] = faster_solverK(V, D, option)
%FN_SOLVER optimizes the bandwidth parameter with kernel function K and normaliztion.
%
%
%
%   See Also: FASTER_SOLVER, FASTER_SOLVER2, FN_SOLVER

%   $ Hyunwoo J. Kim $  $ 2015/11/01 19:17:21 (CST) $

    c1 = option.c1;
    gamma = option.gamma;
    s = option.s0;
    tol = option.tol;
    kernel = option.kernel;
    niter = option.niter;
    
    lambda = diag(D);
    NO = length(lambda); %number of observed variable
    
    ghistory = zeros(niter,1);
    fhistory = zeros(niter,1);
    
    for k=1:niter
        [K,dK] = kernel(s, 1./lambda);
        A = V*diag(K)*V';
        fhistory(k) = lambda'*K - sum(log(lambda.*K))- NO + gamma*norm(A,1)/NO;
        grad = lambda'*dK - sum(dK./K) + gamma/NO*sum(sum(sign(A).*(V*diag(dK)*V')));
        
        grad = grad/NO^2;
        ghistory(k) = grad;

        s_new = s - c1*grad;
        s = s_new;
        
        if s < 0.00
            s = 0.00;
        end    

        fprintf('iter %d, f=%f |grad|=%f, s=%f\n',k, fhistory(k), norm(ghistory(k)), s);
        if norm(grad) <= tol
            ghistory = ghistory(1:k);
            fhistory = fhistory(1:k);
            break;
        end      
    end
    
end

