% Question 4, part b

% This code calculates the external source.
clear all

N = 4;                                          % quadrature order
L = 3;                                          % l in spherical harmonics
source = zeros(21,21);

for l = 0:L                                     
clearvars mu
syms xe eta mu                                  % components of \Omega
stepping = 1;

% compute the integrated terms
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

    f(xe, eta, mu) = term1 * C_lm * P_lm * cos(m * phi);  % Y_{lm}^e
    f2(xe, eta, mu) = term1 * C_lm * P_lm * sin(m * phi);   % Y_{lm}^o

    % perform the integration using LQn quadrature
    [integral_even] = LQnQuadrature(N, f);           % integrate Y_{lm}^e
    [integral_odd] = LQnQuadrature(N, f2);           % integrate Y_{lm}^o
    
    q(stepping) = integral_even;
    s(stepping) = integral_odd;
    stepping = stepping + 1; 
end

i = 2;

% sum over the terms
for m = [0,1:l]
spacing = 0.15;
theta = meshgrid(0:spacing:pi);
mu = cos(theta);
phi = meshgrid(linspace(0, 2*pi, length(theta)));
[Theta, Phi] = meshgrid(0:spacing:pi, linspace(0, 2*pi, length(theta)));

% compute Legendre polynomials derivative term
[derivative] = LegendrePolynomialDerivative(l, m);

% plotting domain for mu
x = mu;

% for Legendre polynomials
prefactor = ((-1).^m) / ((2^l) * factorial(l));
prefactor2 = (1 - x.^2) .^ (m./2);
P_lm = prefactor .* prefactor2 .* eval(derivative);

% for spherical harmonics
prefactor = sqrt((2*l + 1) * factorial(l - m));
prefactor2 = sqrt(4 * pi * factorial(l + m));
Y_lm_e =  term1 .* ((-1).^ m) .* (prefactor ./ prefactor2) .* P_lm .* cos(m .* Phi);
Y_lm_o =  term1 .* ((-1).^ m) .* (prefactor ./ prefactor2) .* P_lm .* sin(m .* Phi);

% perform the summation of the first term
if m ~= 0
    source = source + q(i) * Y_lm_e + s(i) * Y_lm_o;
    i = i + 1;
else
    source = source + q(1) * Y_lm_e;
end
end



end

surf(Theta, Phi, source);
zlabel('Source')
xlabel('Theta (\theta)')
ylabel('Phi (\phi)')
plot_title = sprintf('l_%i', L);
saveas(gcf, plot_title, 'pdf')
