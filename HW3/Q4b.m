% Question 4, part b

% This code calculates the external source.

clear all

N = 6;                                          % quadrature order
l = 8;                                          % l in spherical harmonics

syms xe eta mu                                  % components of \Omega

for m = [0, 1:l]
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

    P_lm = ((-1).^m) / ((2^l) * factorial(l)) .* (1 - mu.^2) .^ (m./2) .* derivative;

    if m == 0
        term1 = sqrt(2-1);
    else
        term1 = sqrt(2);
    end

    theta = acos(mu);
    phi = asin(eta/sin(theta));

    f(xe, eta, mu) = term1 * C_lm * P_lm * cos(m * phi);    % Y_{lm}^e
    f2(xe, eta, mu) = term1 * C_lm * P_lm * sin(m * phi);   % Y_{lm}^o

    % perform the integration using LQn quadrature
    [integral_even] = LQnQuadrature(N, f);           % integrate Y_{lm}^e
    %[integral_odd] = LQnQuadrature(N, f2);           % integrate Y_{lm}^o
    
end

% perform the plotting
spacing = 0.15;
theta = meshgrid(0:spacing:pi);
mu = cos(theta);
phi = meshgrid(linspace(0, 2*pi, length(theta)));
[Theta, Phi] = meshgrid(0:spacing:pi, linspace(0, 2*pi, length(theta)));
