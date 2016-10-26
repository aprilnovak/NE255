% HW4, Question 3
clear all

fontsize = 16;

alpha = 0;
mu = 0.1;
Et = 1.0;
Es = 0.0;
q = 0.0;
h = 0.5;

mesh = 0:h:2;
elem_length = mesh(2);

coordinates = 0:(elem_length / 2):2;

psi(1) = 0; % incoming flux boundary condition

% positive sweep
i = 1;

while i <= (length(coordinates) - 2)
% compute cell-centered value
psi(i + 1) = (q + (2 / (1 + alpha)) * (abs(mu) * psi(i) ./ h)) / (Es + (2 / (1 + alpha)) * abs(mu) / h);

% compute out-going value
psi(i + 2) = (2 / (1 + alpha)) * psi(i + 1) - ((1 - alpha)/(1 + alpha)) * psi(i);

% move to next cell
i = i + 2;
end

% plot the cell-centered values of angular flux
plot(coordinates(2:2:end), psi(2:2:end), '*')
xlabel('Problem Domain', 'FontSize', fontsize)
ylabel('Angular Flux', 'FontSize', fontsize)
