%% Randomly Occurring Deterministic Disturbances Models
%
% Author: Bill Tubbs, February 2022.
%
% Simulation of a class of disturbance models based on concepts and 
% findings from the following paper:
% Title: Duality Between the Control of Processes Subject to Randomly 
%        Occurring Deterministic Disturbances and ARIMA Stochastic 
%        Disturbances
% Authors: J. F. MacGregor, T. J. Harris. J. D. Wright
% Publication: Technometrics, Vol. 26, No. 4, November 1984.
%
clear all; clc

% Random number seed
seed = 0;
rng default

addpath('~/process-observers')
addpath('~/ml-plot-utils/')

results_dir = 'results';
plot_dir = '../images/';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end


%% Generate disturbance signals

% Simulation parameters
Ts = 1; nT = 1000;
t = Ts*(0:nT)';

% Random shock parameters (see MacGregor et al, 1984)
p_shock = 0.01;
sigma_alpha = 1;

% 1. Generate random shock signal
rng(seed)
alpha = sample_random_shocks(size(t), p_shock, sigma_alpha);

% 2. Step Disturbance Process

% Z-transform
Hz = tf([1 0],[1 -1],Ts);
Dt_step = lsim(Hz,alpha,t);

% Difference equation
Hq = tf([1 0],[1 -1],Ts,'Variable','q^-1');
Dt_step_q = lsim(Hq,alpha,t);
assert(max(abs(Dt_step - Dt_step_q)) < 1e-10)

% 3. Ramp Disturbance Process

% Z-transform
Hz = tf([0 0 1],conv([1 -1],[1 -1]),Ts);
Dt_ramp = lsim(Hz,alpha,t);

% Difference equation
Hd = tf([0 0 1],conv([1 -1],[1 -1]),Ts,'Variable','q^-1');
Dt_ramp_d = lsim(Hd,alpha,t);
assert(max(abs(Dt_ramp - Dt_ramp_d)) < 1e-10)

% State-space model
A = [1 1; 0 1]; B = [0; 1]; C = [1 0]; D = 0;
Hss = ss(A,B,C,D,Ts);
Dt_ramp_ss = lsim(Hss,alpha,t);
assert(max(abs(Dt_ramp - Dt_ramp_ss)) < 1e-10)

% 4. Exponential Change Disturbance Process
phi = 0.95;

% Z-transform
Hz = tf([0 0 1],conv([1 -phi],[1 -1]),Ts);
Dt_exp = lsim(Hz,alpha,t);

% Difference equation
Hd = tf([0 0 1],conv([1 -phi],[1 -1]),Ts,'Variable','q^-1');
Dt_exp_d = lsim(Hd,alpha,t);
assert(max(abs(Dt_exp - Dt_exp_d)) < 1e-10)

% State-space model
A = [1 1; 0 phi]; B = [0; 1]; C = [1 0]; D = 0;
Hss = ss(A,B,C,D,Ts);
Dt_exp_ss = lsim(Hss,alpha,t);
assert(max(abs(Dt_exp - Dt_exp_ss)) < 1e-10)


%% Make plot figures

figure(1); clf

ax1 = subplot(4,1,1);
make_tsplot(alpha, t, {'$w_p(k)$'}, [], nan(2), '(a) Random shocks')

ax2 = subplot(4,1,2);
make_tsplot(Dt_step, t, {'$p(k)$'}, [], nan(2), '(b) RODD step disturbance')

ax3 = subplot(4,1,3);
make_tsplot(Dt_ramp, t, {'$p(k)$'}, [], nan(2), '(c) RODD ramp disturbance')

ax4 = subplot(4,1,4);
make_tsplot(Dt_exp, t, {'$p(k)$'}, '$t$', nan(2), ...
    '(d) RODD exponential change disturbance')

linkaxes([ax1 ax2 ax3 ax4], 'x');
set(gcf,'Position',[100 100 400 450])
saveas(gcf,fullfile(plot_dir,'rodd_sim_plots.png'))
save_fig_to_pdf(fullfile(plot_dir,'rodd_sim_plots.pdf'))


%% Combined steps and ramps

rng(seed)
p_shock2 = [p_shock; p_shock];
sigma_alpha2 = [sigma_alpha; 0.01*sigma_alpha];
alpha2 = sample_random_shocks([length(t) 2], p_shock2, sigma_alpha2);

% Z-transform model


% State-space model
A = [1 1; 0 1]; B = [1 0; 0 1]; C = [1 0]; D = 0;
H = ss(A,B,C,D,Ts);
Dt_step_ramp = lsim(H,alpha2,t);

% Plot
figure(2); clf
make_tsplot(Dt_step_ramp, t, {'$p(k)$'})
set(gcf,'Position',[100 625 400 150])
saveas(gcf,fullfile(plot_dir,'rodd_sim_plot2.png'))
save_fig_to_pdf(fullfile(plot_dir,'rodd_sim_plot2.pdf'))


%% Add measurement noise to data

rng(seed)
sigma_M = [0.15 0.15 100 5];
V = sigma_M .* randn(nT+1, 4);

% Add to signals
alpha_M = alpha + V(:, 1);
Dt_step_M = Dt_step + V(:, 2);
Dt_ramp_M = Dt_ramp + V(:, 3);
Dt_exp_M = Dt_exp + V(:, 4);

figure(3); clf

ax1 = subplot(4,1,1);
make_tsplot(alpha_M, t, {'$w_p(k)$'}, [], nan(2), '(a) Random shocks')

ax2 = subplot(4,1,2);
make_tsplot(Dt_step_M, t, {'$p(k)$'}, [], nan(2), '(b) RODD step disturbance')

ax3 = subplot(4,1,3);
make_tsplot(Dt_ramp_M, t, {'$p(k)$'}, [], nan(2), '(c) RODD ramp disturbance')

ax4 = subplot(4,1,4);
make_tsplot(Dt_exp_M, t, {'$p(k)$'}, '$t$', nan(2), ...
    '(d) RODD exponential change disturbance')

linkaxes([ax1 ax2 ax3 ax4], 'x');

% Save figure to image file
set(gcf,'Position',[500 100 400 450])
saveas(gcf,fullfile(plot_dir,'rodd_sim_plots_m.png'))
save_fig_to_pdf(fullfile(plot_dir,'rodd_sim_plots_m.pdf'))

% Save data to file
data = table(t, alpha, Dt_step, Dt_ramp, Dt_exp, alpha_M, Dt_step_M, ...
    Dt_ramp_M, Dt_exp_M);
filename = 'rod_sim_data.csv';
writetable(data, fullfile(results_dir, filename))
fprintf("Simulation data saved to '%s'\n", filename)

