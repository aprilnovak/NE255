% QUESTION 2 - NE 255 hw 2
clear all

N = 8;                                          % quadrature order
syms xe eta mu                                  % components of \Omega

% Parts A and B
%f(xe, eta, mu) = sqrt(xe^2 + eta^2 + mu^2);    % function to integrate

% Part C (uncomment to see the result)
%f(xe, eta, mu) = mu;
%f(xe, eta, mu) = mu^2;
f(xe, eta, mu) = cos(mu);


[wt, valid_full] = LQnQuadrature(N);

% perform the integration over all octants
integral = 0;
j = 1;
for i = 1:(N*(N+2))
    integral = integral +  wt(j) * f(valid_full(i,1), valid_full(i,2), valid_full(i,3));
    if (mod(i, 8) == 0)
        j = j + 1;
    end
end

% multiply result by pi/2 to scale
integral = eval(integral) * pi/2;
disp(sprintf('S-%i quadrature gives an integral of: %.8f', N, integral));
