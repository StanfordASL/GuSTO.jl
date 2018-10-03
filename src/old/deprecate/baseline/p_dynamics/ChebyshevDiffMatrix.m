%   ChebyshevMatrices.m: form cosine-spaced vector of independent variable, 
%                           Chebyshev differentiation and integration matrices
%
%   Ross Allen, ASL, Stanford University
%   Michael Colonno, ADL, Stanford University 
%
%   Started:        6/27/11
%   Last Updated:   Feb 13, 2014
%
%   Inputs:         N               order of polynomial approximation 
%                   t               reverse(original) time points 1 ---> -1 
%
%   Outputs:        D               differentiation matrix (non-dim)
%                   err             nonzero if an error occurs
%
%
%   References: 
%              Q. Gong, I. M. Ross, F. Fahroo, "A Chebyshev Pseudospectral
%               Method for Nonlinear Constrained Optimal Control Problems",
%               Joint 48th IEEE Conference on Decision and Control and 28th
%               Chinese Control Conference, Shanghai, P.R. China, December
%               16-18, 2009
%
%              L. Trfethen, "Spectral Methods in MATLAB", SIAM 2000,
%               Page 128
%
%   Note:
%               Francisco Capistran found problems/inaccuracies in using
%               the integration matrix. Use with caution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = ChebyshevDiffMatrix(N,t)

% differentiation matrix
D = zeros(N+1);

%compute using t (reverse ordered)
for k = 0:N
    i_k = k+1;
    if k == 0 || k==N
        c_k = 2;
    else
        c_k = 1;
    end
    
    for j = 0:N
        i_j = j+1;
        if j == 0 || j==N
            c_j = 2;
        else
            c_j = 1;
        end
    
        if (j == k)     
            D(i_k,i_j) = -0.5*t(i_k)/(1 - t(i_k)^2);       
        else
            D(i_k,i_j) = (c_k/c_j)*((-1)^(j+k))/(t(i_k) - t(i_j));
        end
        
    end
end

% fix corners
D(1,1) = (2*(N^2) +1)/6; D(N+1,N+1) = -D(1,1);

% adjust for forward ordered time
D = -D; 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%