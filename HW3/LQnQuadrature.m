function [wt, valid_full] = LQnQuadrature(N)
% This function provides the weights and quadrature points for a given N
switch N
    case 4
        mu_n = [0.3500212, 0.8688903];
        wt1 = 1/3;
        wt = [wt1, wt1, wt1];
    case 6
        mu_n = [0.2666355, 0.6815076, 0.9261808];
        wt1 = 0.1761263;
        wt2 = 0.1572071;
        wt = [wt1, wt2, wt1, wt2, wt2, wt1];
    case 8
        mu_n = [0.2182179, 0.5773503, 0.7867958, 0.9511897];
        wt1 = 0.1209877;
        wt2 = 0.0907407;
        wt3 = 0.0925926;
        wt = [wt1, wt2, wt2, wt1, wt2, wt3, wt2, wt2, wt2, wt1];
    otherwise
        disp('Non-supported quadrature order.');
end

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

valid = zeros(N*(N+2)/8, 3);
k = 1;
% check to see which quadrature combinations are valid (2-norm equals 1)
for i = 1:total_possible
    norm = sqrt(total(i,1)^2 + total(i,2)^2 + total(i,3)^2);
    if (abs(norm - 1) < 1e-2)
        valid(k,:) = total(i,:);
        k = k + 1;
    end
end

% determine all quadrature points in the other octants
valid_full = zeros(N*(N+2), 3);
i = 1;      
j = 1;
k = 1;
p = 1;

for row = 1:(N*(N+2))
    valid_full(row,:) = [i * valid(p,1), j * valid(p,2), k * valid(p,3)];
    if k == 1
        k = -1;
    else
        k = 1;
        if j == 1;
            j = -1;
        else
            j = 1;
            if i == 1
                i = -1;
            else
                i = 1;
            end
        end
    end
    
    if (mod(row, 8) == 0)
        p = p + 1;
    end
end


end

