% HW4, Question 3
clear all

fontsize = 16;

alpha = 0;
mu = 0.1;
Et = 1.0;
Es = 0.0;
h = 0.1;

mesh = 0:h:2;
elem_length = mesh(2);
coordinates = 0:(elem_length / 2):2;

q = 0.25 .* coordinates;

if (elem_length > 2 * mu / Et)
    disp('Warning: Negative flux!')
end


psi(1) = 2.0; % incoming flux boundary condition

% positive sweep
i = 1;

while i <= (length(coordinates) - 2)
% compute cell-centered value
psi(i + 1) = (q(i+1) + (2 / (1 + alpha)) * (abs(mu) * psi(i) ./ elem_length)) / (Et + (2 / (1 + alpha)) * abs(mu) / elem_length);

% compute out-going value
psi(i + 2) = (2 / (1 + alpha)) * psi(i + 1) - ((1 - alpha)/(1 + alpha)) * psi(i);

% move to next cell
i = i + 2;
end

% negative sweep - store the outgoing flux to apply periodic BC
i = length(coordinates);
psi_n(i) = psi(end);

while i > 2
    
    % compute the cell-centered flux
    psi_n(i - 1) = (q(i-1) + (2/(1 - alpha)) * (mu * psi_n(i) / elem_length)) / (Et + 2 * mu / (elem_length * (1 - alpha)));
    
    % compute the outgoing face flux
    psi_n(i - 2) = (2 / (1 - alpha)) * psi_n(i-1) - ((1+alpha)/(1-alpha)) * psi_n(i);
    
    i = i - 2;
end




% plot the cell-centered values of angular flux
plot(coordinates(2:2:end), psi(2:2:end), '*', coordinates(2:2:end), psi_n(2:2:end), 'r*')
%plot(coordinates(2:2:end), psi_n(2:2:end), 'r*')
xlabel('Problem Domain', 'FontSize', fontsize)
ylabel('Angular Flux', 'FontSize', fontsize)
