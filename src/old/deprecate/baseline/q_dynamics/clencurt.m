% CLENCURT   nodes x (Chebyshev points) and weights w 
%            for Clenshaw-Curtis quadrature
%
%   Reference: L. Trfethen, "Spectral Methods in MATLAB", SIAM 2000,
%              Page 128

function [x,w_net] = clencurt(K)

theta = pi*(0:K)/K;
x = cos(theta); %1 ---> -1 (t)
x = fliplr(x); % -1 ---> 1 (t')

if mod(K,2) == 0
    w = zeros(1,K/2+1); %s:0 ---> N/2
    w(1) = 1/(K^2-1); %w_0
    w_pre = (1/2) + (1/2)*(1/(1-K^2))*cos(pi*(1:K/2)); % s:1--->N/2
    j = (1:K/2-1)';
    for s = 1:K/2
        w(s+1) = (4/K)*w_pre(s) + (4/K)*sum((1./(1-4*(j.^2))).*cos(2*pi*j*s/K));
    end
    w_net = [w, fliplr(w(1:end-1))];
end

w_net = w_net';

end