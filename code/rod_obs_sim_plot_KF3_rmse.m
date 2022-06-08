% Make plots of RMSE for different Kalman Filter tunings from
% simulation results summary file
%
% This is for the Figure 'Tuning of Kalman filter to minimise
% estimation errors' in thesis report.
%

clear all; clc

addpath('~/ml-plot-utils')

% Main folder where all simulation results are saved
results_dir = 'results';

% Simulation group name
sim_name = "rod_obs_sim1_3KF_Q";

% Create main directory for these simulations
base_dir = fullfile(results_dir, sim_name);
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Sub-directory to store plot images
plot_dir = fullfile(base_dir, 'plots');
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end


%% Load simulation summary results from csv file

filename = "rod_obs_sim_outputs_summary.csv";

% Specify column names for which values should be identical
id_cols = {'t_stop', 'Ts', 'nT', 'nu', 'ny', 'n', ...
    'epsilon', 'sigma_M', ...
    'sigma_wp_1', 'sigma_wp_2', ...
    'KF3_R', 'KF3_Q_1_1' ...
};

% Specify column names for which values should be unique
uniq_cols = {'seed', 'KF1_Q_2_2', 'KF2_Q_2_2', 'KF3_Q_2_2'};

%     For 2x2 system
%         filename = "rod_test_obs2_2_4_summary_210712b";
%         
%         % Specify column names for which values should be identical
%         id_cols = {'t_stop', 'Ts', 'nT', 'nu', 'ny', 'n', ...
%             'epsilon_1', 'epsilon_2', 'sigma_M_1', 'sigma_M_2' ...
%             'sigma_w_1_1', 'sigma_w_2_1', 'sigma_w_1_2', 'sigma_w_2_2', ...
%             'KF3_R_1_1', 'KF3_R_2_2', ...
%             'KF3_Q_1_1', 'KF3_Q_2_2' ...
%         };
% 
%         % Specify column names for which values should be unique
%         uniq_cols = {'KF1_Q_3_3', 'KF1_Q_4_4', 'KF2_Q_3_3', 'KF2_Q_4_4', 'KF3_Q_3_3', 'KF3_Q_4_4'};
%         
%         % X-axis variable to plot
%         x_label = 'KF3_Q_3_3';
% 
%         % Choose Y-axes variables to plot
%         y_labels = {'RMSE_x1_est_KF3', 'RMSE_x2_est_KF3', 'RMSE_x3_est_KF3', 'RMSE_x4_est_KF3'};

results_table = readtable(fullfile(base_dir, 'results', filename));
fprintf("Existing results loaded from file: %s\n", filename)

% Choose which results to use
% sigma_M = [0.1; 0.1];
% selected = all(results_table{:, {'sigma_M_1', 'sigma_M_2'}} == repmat([0.1; 0.1]', size(results_table, 1), 1), 2);
% results_table = results_table(selected, :);

[n_rows, n_cols] = size(results_table);
fprintf("Table size: (%d, %d)\n", n_rows, n_cols)

col_names = results_table.Properties.VariableNames;
n_obs = 0;
obs_cols = {};
n_mse = 0;
mse_cols = {};
for i = 1:numel(col_names)
    name = col_names{i};
    if startsWith(name, 'obs_')
        n_obs = n_obs + 1;
        obs_cols{n_obs} = name;
    elseif startsWith(name, 'RMSE_')
        n_mse = n_mse + 1;
        mse_cols{n_mse} = name;
    end
end

% Check all simulations are with same observers and 
% other parameters match or are unique (for each seed)
% (This loop will only execute if there is more than one row)
for i = 2:n_rows

    obs_labels = results_table{i, obs_cols};
    assert(isequal(obs_labels, results_table{1, obs_cols}));

    match_params = results_table{i, id_cols};
    assert(isequal(match_params, results_table{1, id_cols}));

    uniq_params = results_table{i, uniq_cols};
    assert(~isequal(uniq_params, results_table{1, uniq_cols}));
    
end

obs_labels = unique(results_table{:, obs_cols});
fprintf("Observers: %s\n", strjoin(obs_labels));

% Iterate over each set of results
seed_values = unique(results_table.seed);
fprintf("Seeds: %s\n", strjoin(string(seed_values)));


%% Plot RMSEs for one selected simulation

% X-axis variable to plot
x_label = 'KF3_Q_2_2';

% Choose Y-axes variables to plot
y_labels = {'RMSE_y_est_KF3', 'RMSE_x_est_2_KF3'};
y_plot_labels = {'$\mathrm{RMSE}(\hat{Y}(N))$', '$\mathrm{RMSE}(\hat{X}_2(N))$'};
title_strings = { ...
    '(a) Root-mean-squared errors of estimates of the output $y(k)$', ...
    '(b) Root-mean-squared errors of estimates of state $x_2(k)$'};

figure(10); clf

% Plot with one selected seed - this should match the seed used
% to tune the other observers.
seed = 2;

selected_results = results_table(results_table.seed == seed, :);
disp(selected_results)

n = selected_results{1, 'n'};
n_plots = numel(y_labels);
plot_data = selected_results(:, [{x_label} y_labels]);
plot_data = sortrows(plot_data);
[y_mins, ind] = min(plot_data{:, y_labels})
x_argmins = plot_data{ind, x_label}

x_values = plot_data{:, x_label};
for i = 1:n_plots
    subplot(n_plots, 1, i);
    y_values = plot_data{:, y_labels{i}};
    semilogx(x_values, y_values, 'o-', 'Linewidth', 2, 'Markersize', 4);
    ylim([0 inf]);
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$\sigma_{w_p}^2$', 'Interpreter', 'Latex')
    ylabel(y_plot_labels{i}, 'Interpreter', 'Latex')
    grid on
    title(title_strings{i}, 'Interpreter', 'Latex') 
end

% Standard size is [560  420]
set(gcf, 'Position', [100 730 448 336])

filename = sprintf('rod_obs_sim1_plot_KF3_tuning');
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Plot of lines for all seeds simulated

fig = figure(11); clf

axs = cell(1, n_plots);
for i = 1:n_plots
    axs{i} = subplot(n_plots, 1, i);
end

n_seeds = numel(seed_values);
labels = cell(1, n_seeds);
x_argmins = nan(n_seeds, n);
for i_seed = 1:n_seeds

    seed = seed_values(i_seed);
    selected_results = results_table(results_table.seed == seed, :);

    n = selected_results{1, 'n'};
    n_plots = numel(y_labels);
    plot_data = selected_results(:, [{x_label} y_labels]);
    plot_data = sortrows(plot_data);
    [y_mins, ind] = min(plot_data{:, y_labels});
    x_argmins(i_seed, :) = plot_data{ind, x_label}';

    x_values = plot_data{:, x_label};
    for i = 1:n_plots
        subplot(n_plots, 1, i);
        y_values = plot_data{:, y_labels{i}};
        semilogx(x_values, y_values, 'o-', 'Linewidth', 2, 'Markersize', 4);
        hold on
    end
    labels{i_seed} = sprintf("%d%", seed);
end

for i = 1:n_plots
    axes(axs{i})
    ylim([0 inf]);
    hold on
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$\sigma_{w_p}^2$', 'Interpreter', 'Latex')
    ylabel(y_plot_labels{i}, 'Interpreter', 'Latex')
    grid on
    title(title_strings{i}, 'Interpreter', 'Latex')
end

% leg = legend(labels, 'Interpreter', 'Latex');
% leg.Position(1) = 0.9;
% leg.Position(2) = 0.4;

% Standard size is [560  420]
set(fig, 'Position', [100 350 448 336])

filename = sprintf('rod_obs_sim1_plot_KF3_tuning_seeds');
save_fig_to_pdf(fullfile(plot_dir, filename))

fprintf("Best overall values:\n")
disp(median(x_argmins))

return


%% Similar plots for KF1 and KF2 on 2x2 system

% Choose axes to plot
x_label = 'KF1_Q_3_3';
y_labels = {'RMSE_x1_est_KF1', 'RMSE_x2_est_KF1', 'RMSE_x3_est_KF1', 'RMSE_x4_est_KF1'};
n_plots = numel(y_labels);
plot_data = results_table(:, [{x_label} y_labels]);
plot_data = sortrows(plot_data);

figure(11); clf
x_values = plot_data{:, x_label};
for i = 1:n_plots
    subplot(n_plots, 1, i);
    y_values = plot_data{:, y_labels{i}};
    subplot(n,1,i);
    semilogx(x_values, y_values, 'o-')
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$Q_{3,3}$', 'Interpreter', 'Latex')
    ylabel('RMSE', 'Interpreter', 'Latex')
    grid on
    title(sprintf('(%s) Root-mean-squared errors of estimates of $x_%d(k)$', char(96+i), i), ...
    'Interpreter', 'Latex')
end

filename = sprintf('rod_obs_test2_%d_%d_KF1_Qtest.png', p_case, test);
saveas(gcf, fullfile(plot_dir, filename))


% Choose axes to plot
x_label = 'KF2_Q_3_3';
y_labels = {'RMSE_x1_est_KF2', 'RMSE_x2_est_KF2', 'RMSE_x3_est_KF2', 'RMSE_x4_est_KF2'};
n_plots = numel(y_labels);
plot_data = results_table(:, [{x_label} y_labels]);
plot_data = sortrows(plot_data);

figure(12); clf
x_values = plot_data{:, x_label};
for i = 1:n_plots
    subplot(n_plots, 1, i);
    y_values = plot_data{:, y_labels{i}};
    subplot(n,1,i);
    semilogx(x_values, y_values, 'o-')
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$Q_{3,3}$', 'Interpreter', 'Latex')
    ylabel('RMSE', 'Interpreter', 'Latex')
    grid on
    title(sprintf('(%s) Mean-squared errors of estimates of $x_%d(k)$', char(96+i), i), ...
    'Interpreter', 'Latex')
end

filename = sprintf('rod_obs_test2_%d_%d_KF2_Qtest.png', p_case, test);
saveas(gcf, fullfile(plot_dir, filename))

