% HW4, Question 3
clear all

% logistical things
fontsize = 16;

% material data
mu = 0.2;
Et = 1.0;
Es = 0.5;
tolerance = 0.05;
norm = 1; % intial arbitrary value
num_iterations = 0;
max_iterations = 100;

for alpha = [0.0]
    plot_num = 1;

    for h = [0.08]

        mesh = 0:h:2;
        elem_length = mesh(2);
        coordinates = 0:(elem_length / 2):2;

        % initial guess for the scattering source is simply set to the 
        % external source (equivalent to guessing psi = 0.0)
        q = 1.0 .* coordinates;

        % iterate
        while ((norm > tolerance) && (num_iterations < max_iterations))
            
            if (elem_length > 2 * mu / Et)
                sprintf('Warning: Possibly negative flux!')
            end

            % initialize
            psi = zeros(1, length(coordinates));
            psi_n = zeros(1, length(coordinates));
            phi = zeros(1, length(coordinates(2:2:end)));

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
            % plot(coordinates(2:2:end), psi(2:2:end), '*', coordinates(2:2:end), psi_n(2:2:end), 'r*')
            % legend('Flux for Positive \mu', 'Flux for Negative \mu')
            % xlabel('Problem Domain', 'FontSize', fontsize)
            % ylabel('Cell-Centered Angular Flux', 'FontSize', fontsize)
            % saveas(gcf, 'AngularFluxh08', 'jpeg')
            % close all


            % apply "quadrature rule" to get scalar flux
            %phi = psi(2:2:end) + psi_n(2:2:end); % equal weights of 1

            %plot_num = plot_num + 1;
       
            % update the scattering source
            q_new = zeros(1, length(psi));
            q_new = Es .* (psi + psi_n);

            % test for convergence
            norm = 0;
            for u = 1:length(q_new)
                norm = norm + (q_new(u) - q(u))^2;
            end
            norm = sqrt(norm);
            
            q = q_new;
            num_iterations = num_iterations + 1;
            if num_iterations >= max_iterations
                disp('Maximum number of iterations reached!')
            end
gi
            % compute scalar flux
            phi = 0.5 .* (psi(2:2:end) + psi_n(2:2:end));
            
        end
    end
    
end

% plot the scalar flux
% subplot(3, 2, plot_num)
plot(coordinates(2:2:end), phi, '*')
hold on
grid on
legend(sprintf('h = %.2f', h))
xlabel('Problem Domain', 'FontSize', fontsize-2)
ylabel('Scalar Flux', 'FontSize', fontsize-2)
%saveas(gcf, sprintf('ScalarFlux%i', 100*alpha), 'jpeg')
%close all