% QUESTION 2 - NE 255 hw 2
clear all
N = 4;          % quadrature order

switch N
    case 4
        mu_n = [0.3500212, 0.8688903];
    otherwise
        disp('Incorrect quadrature order.');
end
mu_n = [1 2 3];

% permute to find the total number of possible quadrature points
i = 1;              % index for first
j = 1;              % index for second
k = 1;              % index for third
total_possible = length(mu_n)^3;

total = zeros(total_possible, 3);
total(1,:) = [mu_n(1), mu_n(1), mu_n(1)];

for row = 1:total_possible
    total(row,:) = [mu_n(i), mu_n(j), mu_n(k)];
    if (k ~= length(mu_n))
        k = k + 1;
    else
        k = 1;
        if (j ~= length(mu_n))
            j = j + 1;
        else
            j = 1;
            if (i ~= length(mu_n))
                i = i + 1;
            else
                i = 1;
            end
        end
    end
end