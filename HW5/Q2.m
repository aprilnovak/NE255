% HW4, Question 3
clear all

% logistical things
fontsize = 16;

% discrete angles
mu02 = 0.2;
mu05 = 0.5;
mu07 = 0.7;

% group cross sections
Et = [0.5, 0.8, 1.0];
Es = [0.1+0.3+0.1, 0.1+0.3, 0.1+0.3];

% convergence specifications
tolerance = 1e-6;
num_iterations = 0;
max_iterations = 1000;

alpha = 0.5;
h = 0.10;
G = 3;

mesh = 0:h:2;
elem_length = mesh(2);
coordinates = 0:(elem_length / 2):2;

% initial guess for the scattering source for __each__ group. After each 
% within-group calculation, this source will be specialized depending
% on the energy group.

for g = 1:G
        q = 1.0 .* coordinates;
        norm = 1;
        % each within-group calculation finds the angular flux for that group

        % initial guess for the scattering source is simply set to the 
        % external source (equivalent to guessing psi = 0.0)

        % iterate
        while ((norm > tolerance) && (num_iterations < max_iterations))

            % initialize
            psi02 = zeros(1, length(coordinates));
            psi_n02 = zeros(1, length(coordinates));
            psi05 = zeros(1, length(coordinates));
            psi_n05 = zeros(1, length(coordinates));
            psi07 = zeros(1, length(coordinates));
            psi_n07 = zeros(1, length(coordinates));

            psi02(1) = 2.0; % incoming flux boundary condition
            psi05(1) = 2.0; % incoming flux boundary condition
            psi07(1) = 2.0; % incoming flux boundary condition

            % positive sweep
            i = 1;

            while i <= (length(coordinates) - 2)
                % compute cell-centered value
                psi02(i + 1) = (q(i+1) + (2 / (1 + alpha)) * (abs(mu02) * psi02(i) ./ elem_length)) / (Et(g) + (2 / (1 + alpha)) * abs(mu02) / elem_length);
                psi05(i + 1) = (q(i+1) + (2 / (1 + alpha)) * (abs(mu05) * psi05(i) ./ elem_length)) / (Et(g) + (2 / (1 + alpha)) * abs(mu05) / elem_length);
                psi07(i + 1) = (q(i+1) + (2 / (1 + alpha)) * (abs(mu07) * psi07(i) ./ elem_length)) / (Et(g) + (2 / (1 + alpha)) * abs(mu07) / elem_length);

                % compute out-going value
                psi02(i + 2) = (2 / (1 + alpha)) * psi02(i + 1) - ((1 - alpha)/(1 + alpha)) * psi02(i);
                psi05(i + 2) = (2 / (1 + alpha)) * psi05(i + 1) - ((1 - alpha)/(1 + alpha)) * psi05(i);
                psi07(i + 2) = (2 / (1 + alpha)) * psi07(i + 1) - ((1 - alpha)/(1 + alpha)) * psi07(i);
   
                % move to next cell
                i = i + 2;
            end

            % negative sweep - store the outgoing flux to apply periodic BC
            i = length(coordinates);
            psi_n02(i) = psi02(end);
            psi_n05(i) = psi05(end);
            psi_n07(i) = psi07(end);

            while i > 2

                % compute the cell-centered flux
                psi_n02(i - 1) = (q(i-1) + (2/(1 - alpha)) * (mu02 * psi_n02(i) / elem_length)) / (Et(g) + 2 * mu02 / (elem_length * (1 - alpha)));
                psi_n05(i - 1) = (q(i-1) + (2/(1 - alpha)) * (mu05 * psi_n05(i) / elem_length)) / (Et(g) + 2 * mu05 / (elem_length * (1 - alpha)));
                psi_n07(i - 1) = (q(i-1) + (2/(1 - alpha)) * (mu07 * psi_n07(i) / elem_length)) / (Et(g) + 2 * mu07 / (elem_length * (1 - alpha)));

                % compute the outgoing face flux
                psi_n02(i - 2) = (2 / (1 - alpha)) * psi_n02(i-1) - ((1+alpha)/(1-alpha)) * psi_n02(i);
                psi_n05(i - 2) = (2 / (1 - alpha)) * psi_n05(i-1) - ((1+alpha)/(1-alpha)) * psi_n05(i);
                psi_n07(i - 2) = (2 / (1 - alpha)) * psi_n07(i-1) - ((1+alpha)/(1-alpha)) * psi_n07(i);

                i = i - 2;
            end
            
            % update the scattering source, and add the external source
            q_new = zeros(1, length(psi02));
            q_new = Es(g) .* (1/6) .* (psi02 + psi_n02 + psi05 + psi_n05 + psi07 + psi_n07) + 1.0.*coordinates;

            % test for convergence
            norm = 0;
            for u = 1:length(q_new)
                norm = norm + (q_new(u) - q(u)).^2;
            end
            norm = sqrt(norm);
            
            q = q_new;
            
            num_iterations = num_iterations + 1;
            if num_iterations >= max_iterations
                disp('Maximum number of iterations reached!')
            end
            
        end
        
        % compute scalar flux
        phi = zeros(1, length(psi02(2:2:end)));
        phi = (1/6) .* (psi02(2:2:end) + psi_n02(2:2:end) + psi05(2:2:end) + psi_n05(2:2:end) + psi07(2:2:end) + psi_n07(2:2:end));
       
        % save the scalar flux to the appropriate energy group
        group_flux(g,:) = phi;
        
        % save the angular flux as well
        group_flux02(g,:) = psi02(2:2:end);
        %0.5 * (psi02(2:2:end) + psi_n02(2:2:end));
end    
        
        
        
        
        
        
                % plot the scalar flux
        %subplot(3, 2, plot_num)
%         plot(coordinates(2:2:end), phi, '*')
%         hold on
%         grid on
%         legend(sprintf('h = %.2f', h))
%         xlabel('Problem Domain', 'FontSize', fontsize-2)
%         ylabel('Scalar Flux', 'FontSize', fontsize-2)
%saveas(gcf, 'ScalarFlux_partD', 'jpeg')
%close all
