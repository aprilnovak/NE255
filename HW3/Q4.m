% Question 4, HW 3
clear all

l = 4;
m = 0;
theta = 0:0.01:pi;
mu = cos(theta);

% compute Associated LP     l, m, mu
[P_lm] = AssociatedLegendre(l, m, mu);

% compute the spherical harmonics
prefactor = sqrt((2*l + 1) * factorial(l - m));
prefactor2 = sqrt(4 * pi * factorial(l + m));

phi = linspace(0, 2*pi, length(mu));
Y_lm = ((-1).^ m) .* (prefactor ./ prefactor2) .* P_lm .* exp(1j .* m .* phi);

% try using meshgrid
theta = meshgrid(0:0.01:pi);
mu = cos(theta);
phi = meshgrid(linspace(0, 2*pi, length(theta)));


syms x

if ((l + m) == 0)
    derivative = (x.^2 -1) .^ l;
elseif ((l + m) < 0)
    disp('Incorrect combination of l and m');
else
    derivative = diff((x.^2 - 1).^l, x); % first derivative is always computed
    if (l + m) > 1
        for i = 1:(l + m - 1)
            derivative = diff(derivative, x);
        end
    end
end

% plotting domain for mu
x = mu;

% for Legendre polynomials
prefactor = ((-1).^m) / ((2^l) * factorial(l));
prefactor2 = (1 - x.^2) .^ (m./2);
P_lm = prefactor .* prefactor2 .* eval(derivative);

% for spherical harmonics
prefactor = sqrt((2*l + 1) * factorial(l - m));
prefactor2 = sqrt(4 * pi * factorial(l + m));

Z = mu.*phi;
mesh(((-1).^ m) .* (prefactor ./ prefactor2) .* P_lm .* exp(1j .* m .* phi));
%mesh(P_lm)
