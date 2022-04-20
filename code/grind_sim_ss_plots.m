% Make plots of steady-state operating points

clear all

addpath('~/yaml')
addpath('~/ml-plot-utils')

% Sub-directories used
data_dir = 'data';
results_dir = 'results';
plot_dir = '../images';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Load steady-state simulation results
filename = "ss_results_ss_CL_rc.csv";
ss_results = readtable(fullfile(results_dir, filename));

% Input variable names (these must match 'u1, u2, u3' in data
in_vars = {'SAG_Nc', 'SAG_J', 'BASE_ORE_MIX'};

% Output variable names to plot
out_vars = {'TONNAGE', 'SAG_POW', 'SAG_PUMP_SOL', 'SAG_OF_P80'};
n_out = numel(out_vars);

% Use same y axes limits for both plots
y_lims = {[110 135], [2200 2600], [60 63], [100 110]};

% Load normal operating points of process variables
op_data = yaml.loadFile(fullfile(data_dir, 'op_data.yml'));
op_pts = op_data.op_pts;


%% Plots like grind-curves with nominal ore mix factor

mix_factor = op_pts.(in_vars{3});
plot_data = ss_results(ss_results.u3 == mix_factor, [{'u1', 'u2'}, out_vars]);

u3_values = unique(plot_data{:, 'u1'});  % results are sorted
n_u1 = numel(u3_values);

figure(1); clf
in_var_name = in_vars{2};
x_label = sprintf("%s (%s)", op_data.desc.(in_var_name), ...
        op_data.units.(in_var_name));

axs = cell(1, n_out);
for i = 1:n_out
    out_var_name = out_vars{i};
    
    axs{i} = subplot(2, 2, i);
    labels = cell(1, n_u1+1);
    handles = cell(1, n_u1);
    for j = 1:n_u1
        u1 = u3_values(j);
        data = sort(plot_data{plot_data.u1 == u1, {'u2', out_var_name}});
        handles{j} = plot(data(:, 1), data(:, 2), 'o-');
        hold on
        labels{j} = sprintf("speed = %g", u1);
    end
    ylim(y_lims{i})
    
    c = get(handles{2},'Color');  % Use color of second line
    plot(op_pts.(in_var_name), op_pts.(out_var_name), 'o', ...
        'MarkerFaceColor', c)
    labels{end} = 'op. pt.';
    
    xlabel(escape_latex_chars(x_label))
    label = sprintf("%s (%s)", op_data.desc.(out_var_name), ...
        op_data.units.(out_var_name));    ylabel(escape_latex_chars(label))
    grid on
end

linkaxes([axs{1} axs{2} axs{3} axs{4}], 'x')
legend(labels,'Location','best','Interpreter','Latex')

filename = "grind_sim_ss_plot1";
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Plots like grind-curves with different ore mix factors

speed = op_pts.(in_vars{1});
plot_data = ss_results(ss_results.u1 == speed, [{'u3', 'u2'}, out_vars]);

u3_values = unique(plot_data{:, 'u3'});  % results are sorted
n_u3 = numel(u3_values);

figure(2); clf
in_var_name = in_vars{2};
x_label = sprintf("%s (%s)", op_data.desc.(in_var_name), ...
        op_data.units.(in_var_name));

axs = cell(1, n_out);
for i = 1:n_out
    out_var_name = out_vars{i};
    
    axs{i} = subplot(2, 2, i);
    labels = cell(1, n_u3+1);
    handles = cell(1, n_u3);
    for j = 1:n_u3
        u3 = u3_values(j);
        data = sort(plot_data{plot_data.u3 == u3, {'u2', out_var_name}});
        handles{j} = plot(data(:, 1), data(:, 2), 'o-');
        hold on
        labels{j} = sprintf("mix factor = %g", u3);
    end
    ylim(y_lims{i})

    c = get(handles{1},'Color');  % Use color of 1st line
    plot(op_pts.(in_var_name), op_pts.(out_var_name), 'o', ...
        'MarkerFaceColor', c)
    labels{end} = 'op. pt.';
    
    xlabel(escape_latex_chars(x_label))
    label = sprintf("%s (%s)", op_data.desc.(out_var_name), ...
        op_data.units.(out_var_name));    ylabel(escape_latex_chars(label))
    grid on
end

linkaxes([axs{1} axs{2} axs{3} axs{4}], 'x')
legend(labels,'Location','best','Interpreter','Latex')

filename = "grind_sim_ss_plot2";
save_fig_to_pdf(fullfile(plot_dir, filename))
