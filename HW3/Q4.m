% Question 4, HW 3
clear all

l = 1;
m = 0;
mu_mesh = 0.01;

% compute Associated LP     l, m, mu mesh size
[P_lm] = AssociatedLegendre(l, m, mu_mesh);

% compute the spherical harmonics
prefactor = sqrt((2*l + 1) * factorial(l - m));
prefactor2 = sqrt(4 * pi * factorial(l + m));

phi = linspace(0, 2*pi, length(-1:mu_mesh:1));
Y_lm = ((-1).^ m) .* (prefactor ./ prefactor2) .* P_lm .* exp(1j .* m .* phi);