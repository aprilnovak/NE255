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
Es_mat = [0.1, 0.0, 0.0; 0.3, 0.1, 0.1; 0.1, 0.3, 0.3];

% convergence specifications
tolerance = 1e-6;
max_iterations = 1000;
energy_tolerance = 1e-6;
num_energy_iterations = 0;
max_energy_iterations = 2000;

alpha = 0.5;
h = 0.10;
G = 3;

mesh = 0:h:2;
elem_length = mesh(2);
coordinates = 0:(elem_length / 2):2;

% external sources
qe(1,:) = 1.5 .* ones(1, length(coordinates));
qe(2,:) = 0.0 .* ones(1, length(coordinates));
qe(3,:) = 0.2 .* ones(1, length(coordinates));

% initial guesses for scattering sources
q_energy = zeros(3, length(coordinates));

% initial guess for the scattering source for __each__ group. After each 1
% within-group calculation, this source will be specialized depending
% on the energy group
energy_norm = 1; % arbitrary initial value
while ((energy_norm > energy_tolerance) && (num_energy_iterations < max_energy_iterations))
    for g = 1:G
            
            if num_energy_iterations == 0
                q = 1.0 .* coordinates;
            else
                q = q_energy(g,:);
            end
            
            norm = 1;
            % each within-group calculation finds the angular flux for that group

            % initial guess for the scattering source is simply set to the 
            % external source (equivalent to guessing psi = 0.0)

            % iterate
            num_iterations = 0;
            while ((norm > tolerance) && (num_iterations < max_iterations))

                % initialize
                psi02 = zeros(1, length(coordinates));
                psi_n02 = zeros(1, length(coordinates));
                psi05 = zeros(1, length(coordinates));
                psi_n05 = zeros(1, length(coordinates));
                psi07 = zeros(1, length(coordinates));
                psi_n07 = zeros(1, length(coordinates));

                if g == 1
                    psi02(1) = 0.5/3; % incoming flux boundary condition
                    psi05(1) = 0.5/3; % incoming flux boundary condition
                    psi07(1) = 0.5/3; % incoming flux boundary condition
                end

                % positive sweep
                i = 1;
                while i <= (length(coordinates) - 2)
                    % compute cell-centered value
                    psi02(i + 1) = (q(i+1) + qe(g, i+1) + (2 / (1 + alpha)) * (abs(mu02) * psi02(i) ./ elem_length)) / (Et(g) + (2 / (1 + alpha)) * abs(mu02) / elem_length);
                    psi05(i + 1) = (q(i+1) + qe(g, i+1) + (2 / (1 + alpha)) * (abs(mu05) * psi05(i) ./ elem_length)) / (Et(g) + (2 / (1 + alpha)) * abs(mu05) / elem_length);
                    psi07(i + 1) = (q(i+1) + qe(g, i+1) + (2 / (1 + alpha)) * (abs(mu07) * psi07(i) ./ elem_length)) / (Et(g) + (2 / (1 + alpha)) * abs(mu07) / elem_length);

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
                    psi_n02(i - 1) = (q(i-1) + qe(g, i-1) + (2/(1 - alpha)) * (mu02 * psi_n02(i) / elem_length)) / (Et(g) + 2 * mu02 / (elem_length * (1 - alpha)));
                    psi_n05(i - 1) = (q(i-1) + qe(g, i-1) + (2/(1 - alpha)) * (mu05 * psi_n05(i) / elem_length)) / (Et(g) + 2 * mu05 / (elem_length * (1 - alpha)));
                    psi_n07(i - 1) = (q(i-1) + qe(g, i-1) + (2/(1 - alpha)) * (mu07 * psi_n07(i) / elem_length)) / (Et(g) + 2 * mu07 / (elem_length * (1 - alpha)));

                    % compute the outgoing face flux
                    psi_n02(i - 2) = (2 / (1 - alpha)) * psi_n02(i-1) - ((1+alpha)/(1-alpha)) * psi_n02(i);
                    psi_n05(i - 2) = (2 / (1 - alpha)) * psi_n05(i-1) - ((1+alpha)/(1-alpha)) * psi_n05(i);
                    psi_n07(i - 2) = (2 / (1 - alpha)) * psi_n07(i-1) - ((1+alpha)/(1-alpha)) * psi_n07(i);

                    i = i - 2;
                end

                % update the within-group scattering source
                q_new = zeros(1, length(psi02));
                q_new = Es(g) .* (1/6) .* (psi02 + psi_n02 + psi05 + psi_n05 + psi07 + psi_n07) + q_energy(g,:)/3;

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

            % compute scalar flux (cell-centered)
            phi = zeros(1, length(psi02(2:2:end)));
            phi = (1/6) .* (psi02(2:2:end) + psi_n02(2:2:end) + psi05(2:2:end) + psi_n05(2:2:end) + psi07(2:2:end) + psi_n07(2:2:end));
            
            phi_total = zeros(1, length(psi02));
            phi_total = (1/6) .* (psi02 + psi_n02 + psi05 + psi_n05 + psi07 + psi_n07);
            
            % save the scalar flux to the appropriate energy group
            group_flux(g,:) = phi;
            group_flux_total(g,:) = phi_total;

    end    

    % after updating all of the groups, update the group sources using the
    % scattering matrix.
    for gg = 1:G
        q_energy(gg,:) = Es_mat(gg,1).*group_flux_total(1,:) + Es_mat(gg,2).*group_flux_total(2,:) + Es_mat(gg,3).*group_flux_total(3,:);
    end
    
    
    % compute whether or not the energy iteration has converged
    energy_norm = 0;
    for u = 1:length(q_energy(1,:))
        energy_norm = energy_norm + (q_energy(1,u) - q(1,u)).^2;
    end
    energy_norm = sqrt(energy_norm);
    
    num_energy_iterations = num_energy_iterations + 1;
    if num_energy_iterations >= max_energy_iterations
        disp('Maximum number of energy iterations reached!')
    end

end

        
        
        
        
% plot the scalar flux
plot(coordinates(2:2:end), group_flux(1,:), '*', coordinates(2:2:end), group_flux(2,:), '*', coordinates(2:2:end), group_flux(3,:), '*')
hold on
grid on
legend('Group 1', 'Group 2', 'Group 3')
xlabel('Problem Domain', 'FontSize', fontsize-2)
ylabel('Scalar Flux', 'FontSize', fontsize-2)
saveas(gcf, 'ScalarFlux_3Group', 'jpeg')
