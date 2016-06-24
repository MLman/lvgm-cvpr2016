function [theta, theta_GT]= random_sparse_network(NO, NC, N, density)
%RANDOM_SPARSE_NETWORK generates undirected graph (sparsity pattern of
%precision matrix).
%
%   NO : number of image derived features.
%   NC : number of meta data features (covariates).
%   N  : Total number of varibles including unobserved (latent) variables.
%   N-NO-NC : number of unobserved (latent) variables.
%   theta : observed precision matrix
%   theta_GT : ground truth precision matrix.
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/06/24 12:48:29 (CDT) $

A = sprand( NO, NO, density ); % generate adjacency matrix at random
A = tril( A, -1 );    
A(A<0.5) = 0;
A = A + A';
ttheta = zeros(N);
ttheta(1:NO, 1:NO) = A;
for i=1:NO
    for j=NO+1:NO+NC        
        e = 0.5;
        if i ~= j
            if i <= NO/NC
                ttheta(i,NO+1) = e;
                ttheta(NO+1,i) = e;            
            elseif i <= 2*NO/NC
                ttheta(i,NO+2) = e;
                ttheta(NO+2,i) = e;    
            elseif i <= 3*NO/NC
                ttheta(i,NO+3) = e;
                ttheta(NO+3,i) = e;    
            elseif i <= 4*NO/NC
                ttheta(i,NO+4) = e;
                ttheta(NO+4,i) = e;    
            elseif i <= 5*NO/NC
                ttheta(i,NO+5) = e;
                ttheta(NO+5,i) = e;              
            elseif i <= 6*NO/NC
                ttheta(i,NO+6) = e;
                ttheta(NO+6,i) = e;              
            elseif i <= 7*NO/NC
                ttheta(i,NO+7) = e;
                ttheta(NO+7,i) = e;              
            elseif i <= 8*NO/NC
                ttheta(i,NO+8) = e;
                ttheta(NO+8,i) = e;              
            elseif i <= 9*NO/NC
                ttheta(i,NO+9) = e;
                ttheta(NO+9,i) = e;              
            elseif i <= 10*NO/NC
                ttheta(i,NO+10) = e;
                ttheta(NO+10,i) = e;              
            end
        end
    end   
end
% Ground Truth theta
theta_GT = ttheta+eye(size(ttheta,1));

for i=1:N    
    for j=NO+NC+1:N
        e = .2*abs(randn);    
        if i ~= j
            ttheta(i,j) = e;
            ttheta(j,i) = e;            
        end
    end   
end
theta = make_spd(ttheta);

