close all; clc; clearvars;

% Constants
coeff = 12.27;
tube_diam = 0.0135; % in meters
diam_err = 0.0003; % in meters
V_err = 50; % in Volts

% Reading data from files and converting to appropriate units
files = dir('X-Ray\results.xlsx');
file_name = fullfile(files.name);
diff_data = table2array(readtable(file_name, "VariableNamingRule", "preserve"));
V = diff_data(:,1)*1000; % converting to Volts
lambda = diff_data(:,7);
R1 = diff_data(:,8)/1000; % converting to meters
R2 = diff_data(:,9)/1000; % converting to meters
V_recip = 1./sqrt(V);
V_recip_err = V_err./(V.^(3/2));

% Initial point guess
startPoints = [1,1];

% linear function model
lin_model = 'a*x + b';

% Fitting and goodness of fit
[f1,gof1] = fit(R1, V_recip, lin_model, 'Start', startPoints);
[f2,gof2] = fit(R2, V_recip, lin_model, 'Start', startPoints);
gof1
gof2


% Calculating the interplanar distance of the 1st two planes
d1 = coeff*tube_diam/f1.a
d2 = coeff*tube_diam/f2.a

err_d1 = coeff*(diam_err/f1.a + (tube_diam*gof1.rmse)/(f1.a)^2)
err_d2 = coeff*(diam_err/f2.a + (tube_diam*gof1.rmse)/(f2.a)^2)

residuals1 = V_recip - (f1.a*R1 + f1.b);
chi_sqr1 = sum(((residuals1./V_recip_err).^2));
df1 = length(V_recip) - 2;
chi_sqr_red1 = chi_sqr1/df1;
P_value1 = 1 - chi2cdf(chi_sqr_red1,df1)

residuals2 = V_recip - (f2.a*R2 + f2.b);
chi_sqr2 = sum(((residuals2./V_recip_err).^2));
df2 = length(V_recip) - 2;
chi_sqr_red2 = chi_sqr2/df2;
P_value2 = 1 - chi2cdf(chi_sqr_red2,df1)

% Plotting results
figure (1)
hold on
grid on
plot(f1, R1, V_recip)
xlabel('r_1 (m)')
ylabel('V^{-1/2}')
errorbar(R1, V_recip, V_recip_err, 'vertical','.')
legend('Data',['a = ', num2str(f1.a), '\pm', num2str(gof1.rmse)])

figure (2)
hold on
grid on
plot(f2, R2, V_recip)
xlabel('r_2 (m)')
ylabel('V^{ -1/2}')
errorbar(R2, V_recip, V_recip_err, 'vertical','.')
legend('Data',['a = ', num2str(f2.a), '\pm', num2str(gof2.rmse)])
