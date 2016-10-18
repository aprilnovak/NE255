% QUESTION 2 - NE 255 hw 2
clear all

N = 8;                                          % quadrature order
syms xe eta mu                                  % components of \Omega

% Parts A and B
%f(xe, eta, mu) = sqrt(xe^2 + eta^2 + mu^2);    % function to integrate

% Part C (uncomment to see the result)
%f(xe, eta, mu) = mu;
%f(xe, eta, mu) = mu^2;
f(xe, eta, mu) = cos(mu);

% perform the integration using LQn quadrature
[integral] = LQnQuadrature(N, f);
