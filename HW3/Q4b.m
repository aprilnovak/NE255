% Question 4, part b

% This code calculates the external source, and hence includes the 
% entirety of the code developed for Question 2.

clear all

N = 8;                                          % quadrature order
syms xe eta mu                                  % components of \Omega

% new stuff for computing the odd/even spherical harmonics
l = 3;
m = -1;
C_lm = (-1)^m * sqrt((2*l+1) * factorial(l-m) / (4 * pi * factorial(l+m)));

% compute Associated Legendre polynomials derivative term
if ((l + m) == 0)
    derivative = (mu.^2 -1) .^ l;
elseif ((l + m) < 0)
    disp('Incorrect combination of l and m');
else
    derivative = diff((mu.^2 - 1).^l, mu);
    if (l + m) > 1
        for i = 1:(l + m - 1)
            derivative = diff(derivative, mu);
        end
    end
end

prefactor = ((-1).^m) / ((2^l) * factorial(l));
prefactor2 = (1 - mu.^2) .^ (m./2);
P_lm = prefactor .* prefactor2 .* derivative;

if m == 0
    term1 = sqrt(2-1);
else
    term1 = sqrt(2);
end

theta = acos(mu);
phi = asin(eta/sin(theta));

f(xe, eta, mu) = term1 * C_lm * P_lm * cos(m * phi);    % integration of Y_{lm}^e
f2(xe, eta, mu) = term1 * C_lm * P_lm * sin(m * phi);   % integration of Y_{lm}^o

% perform the integration using LQn quadrature
[integral] = LQnQuadrature(N, f);
