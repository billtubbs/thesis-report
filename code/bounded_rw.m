%% Bounded Random Walk Disturbance Model
%
% Author: Bill Tubbs, March 2022.
%
% Simulation of the Bounded Random Walk stochastic process proposed
% by J. Nicolau:
% Title: Stationary Processes That Look Like Random Walks - The Bounded
%        Random Walk Process in Discrete and Continuous Time
% Author: J. Nicolau
% Publication: Econometric Theory, 18, 2002, 99-118.
%

clear
rng(0)

addpath('~/ml-plot-utils/')
plot_dir = '../images';

% Model (integrator)
A = 1;
B = 1;
C = 1;

% Model parameters
sd_e = 0.4;    % Noise std. dev.
beta = -15;   % k parameter in Nicolau's paper
alpha1 = 3;
alpha2 = 3;
tau = 100;   % initial value x(0)

% Simulation timesteps
N = 1000;

% Set initial state
x = tau;
x_b = x;
p = zeros(1,N);
p_bound = zeros(1,N);

% Noise input
e = randn(1, N);

for i = 1:N

    % Stochastic input variable
    alpha = sd_e * e(i);

    % Unbounded process output
    p(i) = C*x;
    x = A*x + B*alpha;

    % Bounded process output
    p_bound(i) = C*x_b;
    x_b = A*x_b + B*alpha + ...
        brw_reversion_bias(x_b, alpha1, alpha2, beta, tau);

    %fprintf("%6.1f, %6.1f\n", x, x_b)

end

% Simulation

figure(1); clf
t = 0:N-1;
plot(t, p, t, p_bound, 'Linewidth', 2);
ylim(axes_limits_with_margin([p' p_bound'], 0.05));
yline(95, 'k--');
yline(105, 'k--');
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$p(k)$', 'Interpreter', 'latex')
grid on
legend('RW', 'BRW', 'Interpreter', 'latex', 'location', 'best')
set(gcf,'Position',[100 100 400 150])
save2pdf(fullfile(plot_dir, 'brw_sim.pdf'))

% figure(2); clf
% edges = linspace(80, 120, 41);
% histogram(p, edges, 'Normalization', 'pdf'); hold on
% histogram(p_bound, edges, 'Normalization', 'pdf')
% ylabel("Probability density", 'Interpreter', 'latex')
% grid on
% legend('RW', 'BRW')

% % Analyse differences
% dp = diff(p_bound);
% edges = linspace(96, 104, 41);
% p_bins = discretize(p_bound(1:end-1), edges);
% dp_bin_means = nan(1, length(edges));
% for i = 1:length(edges)
%     dp_bin_means(i) = mean(dp(p_bins == i));
% end
% 
% figure(3); clf
% s = scatter(p_bound(1:end-1), dp, 'k.'); hold on
% set(s, 'MarkerEdgeAlpha', 0.25);
% plot(edges, dp_bin_means, 'o-', 'Linewidth', 2)
% set(gca, 'TickLabelInterpreter', 'latex')
% xlabel('p(k)', 'Interpreter', 'latex')
% ylabel('p(k) - p(k-1)', 'Interpreter', 'latex')
% grid on

%% Plot of Bounded Random Walk bias function

figure(4); clf
x = 94:0.1:106;
% a = brw_reversion_bias(x, 0, alpha2, beta, tau);
% plot(x, a, 'Linewidth', 2); hold on
% a = brw_reversion_bias(x, alpha1, 0, beta, tau);
% plot(x, a, 'Linewidth', 2);
a = brw_reversion_bias(x, alpha1, alpha2, beta, tau);
plot(x, a, 'Linewidth', 2);
set(gca, 'TickLabelInterpreter', 'latex')
grid on
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$a(x)$", 'Interpreter', 'latex')
%title("BRW Reversion Function")
set(gcf,'Position',[100 300 300 150])
save2pdf(fullfile(plot_dir, 'brw_a.pdf'))


%% Stationary pdf

dt = 0.02;
x = 94:dt:106;
brw_pdf = brwpdf(x, alpha1, alpha2, beta, tau, sd_e);

% Normalize the pdf
brw_pdf = brw_pdf ./ (sum(brw_pdf) * dt);

figure(5); clf
plot(x, brw_pdf, 'Linewidth', 2)
set(gca, 'TickLabelInterpreter', 'latex')
grid on
xlabel("$p(k)$", 'Interpreter', 'latex')
ylabel("$\mathrm{Pr}(p(k))$", 'Interpreter', 'latex')
ylim(axes_limits_with_margin(brw_pdf));
set(gcf,'Position',[100 525 300 150])
%save2pdf(fullfile(plot_dir, 'brw_pdf.pdf'))


%% Comparison with RW pdf

N = 20;
rw_pdf = normpdf(x, tau, N(1)*sd_e^2);

figure(6); clf
plot(x, rw_pdf, 'Linewidth', 2); hold on
plot(x, brw_pdf, 'Linewidth', 2)
set(gca, 'TickLabelInterpreter', 'latex')
grid on
xlabel("$p(k)$", 'Interpreter', 'latex')
ylabel("$\mathrm{Pr}(p(k))$", 'Interpreter', 'latex')
ylim(axes_limits_with_margin([rw_pdf' brw_pdf']));
labels = {sprintf('RW ($k = %d$)', N), 'BRW ($k = \infty$)'};
legend(labels, 'Interpreter', 'latex', 'location', 'south')
set(gcf,'Position',[100 750 300 150])
save2pdf(fullfile(plot_dir, 'brw_pdf.pdf'))


function a = brw_reversion_bias(x, alpha1, alpha2, beta, tau)
% This is the function 'a(x)' from Nicolau (2002) used in the
% difference equation of the bounded random walk (BRW) 
% (see Eq. 1 in the paper).
% Note: beta is the 'k' parameter in Nicolau's paper.
%
    a = exp(beta) * (exp(-alpha1 * (x - tau)) - exp(alpha2 * (x - tau)));
end


function p = brwpdf(x, alpha1, alpha2, beta, tau, sd_e)
% This is the function 'a(x)' from Nicolau (2002) used in the
% difference equation of the bounded random walk (BRW) 
% (see Eq. 1 in the paper).
% Note: beta is the 'k' parameter in Nicolau's paper.
%
    p = sd_e^(-2) * exp( ...
        -2 * exp(beta) / sd_e^2 * (exp(-alpha1*(x - tau)) / alpha1 + exp(alpha2*(x - tau)) / alpha2) ...
    );
end