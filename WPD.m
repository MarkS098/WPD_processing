close all; clc; clearvars;
% Constants
scale_factor = 59.33;
offset = 68.8;
rho_0 = 5.65; % micro Ohm * cm
R_room = 1.58; % Ohm
R_0 = 1.2; % Ohm
R_w = R_room - R_0;
I_err = 0.005;
V_err = 0.1;
affix_str = 'V.csv';

files = dir(['Measurements\*', affix_str]);
file_names = string(fullfile({files.name}));
temp_name_arr = erase(file_names, affix_str);
temp_name_arr = sort(double(temp_name_arr), 'ascend');
file_names = strcat(string(temp_name_arr), affix_str);

V = temp_name_arr;
I = [0.48, 0.51, 0.54, 0.58, 0.61, 0.64, 0.68]; % 5V, 6V, 7V, 8V, 9V, 10V, 11V
R = V./I;

% Error calculations for rho and R
R_err = R.*sqrt((V_err./V).^2 + (I_err./I).^2);
rho = rho_0*(R./R_0); % micro Ohms
rho_err = rho.*(R_err./R).^2


% Coefficients for wavelength calculation
A = -49852133;
B = 86092018.9;
C = -29983328.35;
D = -14354236.56;
E = 835425.05;
F = 5647432.02;
G = 1863438.86;
H = -2719226.18;
I = 574967.82;

% Wavelength formula according to instructions
lambda = @(n) 3000./sqrt(A + B*n + C*n.^2 + D*n.^3 + E*n.^4 + F*n.^5 + G*n.^6 + H*n.^7 + I*n.^8);

% Refractive index fomrula
n = @(th) sqrt((2/sqrt(3) * sin(th) + 1/2).^2 + 3/4);

% Temperature function and error
T = @(rho) 103 + 38.1.*rho - 0.095.*rho.^2 + (2.48e-4).*rho.^3;
T_k = T(rho)
T_err = T_k.*(38.1*rho_err./rho - 0.19*rho_err./rho + 7.44*rho_err./rho)

% Wein model
lin_model = 'a*x';

% Initial point guess
startPoints = 1;

for i = 1:width(file_names)
    V_data = table2array(readtable(file_names(i), "VariableNamingRule", "preserve"));
    voltage_str = erase(file_names(i), '.csv');
    
    int_data = (V_data(:, 1:2:end));
    angle_data = (V_data(:, 2:2:end));

    max_int = max(int_data, [], 1);
    for j = 1:width(angle_data)
        n_cell{j} = wpdfilt(1.697, 10, n(wpdfilt(0, 20, (offset - angle_data(:,j))/scale_factor)));

        % Determine the maximum length of n_max across all columns
        max_length = max(cellfun(@length, n_cell));

        % Pad each cell of n_max with -1 to match max_length
        n_cell{j}(end+1:max_length) = -1;
    end

    % Convert n_max cell array to a matrix (each cell becomes a column)
    n_matrix = NaN(max_length, length(n_cell));

    for j = 1:length(n_cell)
        n_matrix(1:length(n_cell{j}), j) = n_cell{j}';
    end

    n_matrix(n_matrix == -1) = NaN;
    lambda_init = lambda(n_matrix);
    lambda_max = rmmissing(lambda_init(int_data == max_int));
    lambda_mean_max(i) = mean(lambda_max);

    figure (i)
    hold on
    grid on

    title(['Measurement at: ', voltage_str])
    subtitle(['mean \lambda_{max} = ', num2str(lambda_mean_max(i)), ' (nm)'])
    plot(angle_data, int_data)    
    xlabel('\theta')
    ylabel('Relative intensity (%)')
    xlim([0, max(angle_data, [], 'all')])
    legend
end

% Performing the fits according to a linear model
f1 = fit((1./T_k)', lambda_mean_max' * 1e-9, lin_model, 'Start', startPoints);

figure (i + 1)
plot(f1, 1./T_k, lambda_mean_max * 1e-9)
title(['Weins constant: ' ,num2str(f1.a), ' m\cdotK'])
xlabel('1/T (K^{-1})')
ylabel('\lambda (m)')