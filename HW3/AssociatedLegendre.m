function [P_lm] = AssociatedLegendre(l, m, mu)
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

prefactor = ((-1).^m) / ((2^l) * factorial(l));
prefactor2 = (1 - x.^2) .^ (m./2);
P_lm = prefactor .* prefactor2 .* eval(derivative);

plot(x, P_lm)