% Question 4, HW 3
clear all

l = 1;
m = 1;

theta = meshgrid(0:0.01:pi);
mu = cos(theta);
phi = meshgrid(linspace(0, 2*pi, length(theta)));
[Theta, Phi] = meshgrid(0:0.01:pi, linspace(0, 2*pi, length(theta)));

% compute Legendre polynomials
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
Y_lm = ((-1).^ m) .* (prefactor ./ prefactor2) .* P_lm .* exp(1i .* m .* Phi);

% plot the real component
% mesh(Theta, Phi, real(Y_lm));
% zlabel('Real Component of Ylm')
% xlabel('\theta')
% ylabel('\phi')

% plot the imaginary component
surf(Theta, Phi, imag(Y_lm));
%
