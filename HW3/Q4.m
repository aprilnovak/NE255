% Question 4, HW 3
clear all

j = 1;              % subplot number
real_plot = 0;      % 1 = plot real components, 0 = plot imaginary

for l = [0, 1, 2];
    for m = -l:l
        spacing = 0.15;
        theta = meshgrid(0:spacing:pi);
        mu = cos(theta);
        phi = meshgrid(linspace(0, 2*pi, length(theta)));
        [Theta, Phi] = meshgrid(0:spacing:pi, linspace(0, 2*pi, length(theta)));

        % compute Legendre polynomials derivative term
        [derivative] = LegendrePolynomialDerivative(l, m);

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

        if (real_plot == 1)
            % plot the real component
            subplot(3, 3, j)
            surf(Theta, Phi, real(Y_lm));
            title(sprintf('l = %i, m = %i', l, m))
            zlabel('Re(Ylm)')
            xlabel('Theta (\theta)')
            ylabel('Phi (\phi)')
            plot_title = 'Real_Ylm';
        else
            % plot the imaginary component
            subplot(3, 3, j)
            surf(Theta, Phi, imag(Y_lm));
            title(sprintf('l = %i, m = %i', l, m))
            zlabel('Im(Ylm)')
            xlabel('Theta (\theta)')
            ylabel('Phi (\phi)')
            plot_title = 'Imag_Ylm';
        end

        j = j + 1;
    end
end

saveas(gcf, plot_title, 'pdf')