function [derivative] = LegendrePolynomialDerivative(l, m)

% This function uses a loop to compute the derivative that appears in the
% definition of the Associated Legendre polynomials.

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


end

