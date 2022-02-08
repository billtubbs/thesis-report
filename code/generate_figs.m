%% Produce figures for thesis report

clear all

addpath("../process-observers/")

% Sub-directories used
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


%% PDF of shocks - MacGregor report

% Probability density functions
alpha = linspace(-3, 3, 601);
p = normpdf(alpha);
p0 = p(alpha == 0);

% Figure
figure(1); clf
plot(alpha, p, 'k-', 'Linewidth', 2); hold on
plot([0 0], [p0 p0+5], 'k-', 'Linewidth', 2);
ylim([0 p0+5])
xticks([0])
yticks([0])
xlabel('$\alpha(k)$', 'Interpreter', 'Latex')
ylabel('$P(\alpha(k))$', 'Interpreter', 'Latex')
set(gcf,'Position',[100 100 300 150])
saveas(gcf,fullfile(plot_dir,'alpha-pdf-1.png'))

%% PDF of shocks - Robertson report

% System parameters
b = 100;
sigma1 = 0.01;
sigma_w = [sigma1; sigma1*b];
epsilon = 0.01;
p_gamma = [1-epsilon; epsilon];

% Conditional probability of w_i(k)
wp = linspace(-2, 2, 401);
p1 = prob_w_given_gamma(wp, 0, sigma_w);
p2 = prob_w_given_gamma(wp, 1, sigma_w);

% Figure
figure(2); clf
plot(wp, p1, 'Linewidth', 2); hold on
plot(wp, p2, 'Linewidth', 2);
grid on
xlabel('$w_p(k)$', 'Interpreter', 'Latex')
ylabel('$p(w_p(k)|\gamma(k))$', 'Interpreter', 'Latex')
legend('$\gamma(k)=0$','$\gamma(k)=1$', 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','\sigma_w',sigma1);
text(-1.3, 5.6, s, 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','b\sigma_w',sigma1*b);
text(-1.3, 5, s, 'Interpreter', 'Latex')
set(gcf,'Position',[100 300 300 200])
saveas(gcf,fullfile(plot_dir,'alpha-pdfs-2.png'))

% Figure - multi-variable
figure(3); clf
plot(wp, p1, 'Linewidth', 2); hold on
plot(wp, p2, 'Linewidth', 2);
grid on
xlabel('$w_{p,i}(k)$', 'Interpreter', 'Latex')
ylabel('$p(w_{p,i}(k)|\gamma_i(k))$', 'Interpreter', 'Latex')
legend('$\gamma_i(k)=0$','$\gamma_i(k)=1$', 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','\sigma_i',sigma1);
text(-1.3, 5.6, s, 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','b_i\sigma_i',sigma1*b);
text(-1.3, 5, s, 'Interpreter', 'Latex')
set(gcf,'Position',[100 500 300 200])


% Total probability of w_i(k)
p = prob_w(wp, p_gamma, sigma_w);

figure(4); clf
plot(wp, p, 'Linewidth', 2)
grid on
%ylim([0 45])
xlabel('$w_p(k)$', 'Interpreter', 'Latex')
ylabel('$\Pr(w_p(k))$', 'Interpreter', 'Latex')
text_array = { ...
    sprintf('$$%s=%g$$','\sigma_{w_p}', sigma1), ...
    sprintf('$$%s=%g$$', 'b', b), ...
    sprintf('$%s=%g$','\epsilon', epsilon(1)) ...
};
text(-1.5, 25, text_array);
set(gcf,'Position',[400 100 300 200])
saveas(gcf,fullfile(plot_dir,'alpha-pdf-4.png'))
saveas(gcf,fullfile(plot_dir,'alpha-pdf-4.svg'))
saveas(gcf,fullfile(plot_dir,'alpha-pdf-4.eps'))

% Calculate how much density
wp = linspace(-10, 10, 1001);  % use a finer scale
p = prob_w(wp, p_gamma, sigma_w);
i_center = find(wp == 0);

% Conditional probability of gamma given w_i(k)
dist1 = makedist('Normal','mu',0,'sigma',0.01);
dist2 = makedist('Normal','mu',0,'sigma',1);
dx = 0.04;
cum_density_total(dx, dist1, dist2, p_gamma);
fun = @(dx) (cum_density_total(dx, dist1, dist2, p_gamma) - 0.99)^2;
fun(dx)
dx = fminsearch(fun, dx)

w = linspace(-3, 3, 601);
p1 = prob_gamma_given_w(0, w, p_gamma, sigma_w);
p2 = prob_gamma_given_w(1, w, p_gamma, sigma_w);

figure(5); clf
plot(w,p1,'Linewidth',2); hold on
plot(w,p2,'Linewidth',2); grid on
xlabel('$w_{p,i}(k)$', 'Interpreter', 'Latex')
ylabel('$p(\gamma_i(k)|w_{p,i}(k))$', 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','\sigma_i', sigma1);
text(-2.4, 0.62, s, 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','b_i\sigma_i', sigma1*b);
text(-2.4, 0.54, s, 'Interpreter', 'Latex')
s = sprintf('$%s=%g$','\epsilon_i', p_gamma(1));
text(-2.4, 0.46, s, 'Interpreter', 'Latex')
legend('$p(\gamma_i(k)=0|w_{p,i}(k))$', ...
    '$p(\gamma_i(k)=1|w_{p,i}(k))$', 'Interpreter', 'Latex')
set(gcf,'Position',[400 300 400 200])


% %% Binomial distribution
% 
% d = [5 10 20 40];
% d_max = 5;
% p = 0.01;
% Pr = nan(numel(d), d_max+1);
% for i=1:numel(d)
%     Pr(i, 1:d_max+1) = binopdf(0:d_max, d(i), p);
% end
% Pr
% 
% figure(6); clf
% labels = {};
% for i=1:numel(d)
%     plot(0:d_max, Pr(i, :), 'o-'); hold on
%     labels{i} = sprintf("d = %d", d(i));
% end
% xlabel('Number of shocks')
% ylabel('Probability')
% grid on
% legend(labels)
% set(gcf,'Position',[400 500 300 200])
% 
% 
% %% PDF of prob_y_given_gamma
% 
% Ts = 0.5;
% 
% % System parameters
% epsilon = 0.01;
% epsilon = [epsilon 1-epsilon];
% sigma_w = [1; 0.01];
% 
% % Load system matrices
% sys_rodin_step
% 
% % Initial condition
% x0 = zeros(n, 1);
% 
% k = 5;
% [F, G] = pred_mats(A,B,C,D,k,k);
% Fk = F(end,:);
% Gk = G(end,:);
% 
% t = Ts*(0:k-1)';
% U = zeros(k, 1);
% U(t >= 5) = 1;
% Wp_exp = zeros(k, 1);
% U_combined = reshape([U Wp_exp]', k*2, 1);
% 
% % Get all combinations of shock sequences
% combs = combinations_lte(k, 10);
% 
% % Expected value of y(k) (same for all sequences)
% y_exp = Fk*x0 + Gk*U_combined;
% 
% % Split Gk into two to separate Uu and Uwp
% Gku = Gk(1:2:end-1);
% Gkwp = Gk(2:2:end);
% 
% % Covariances of y(k) for each sequence Gamma(k)
% 
% % Probability of sequences Gamma(k)
% p_seqs = prob_Gammas(S, epsilon);
% 
% return
% 
% %y_cov = 
% 
% sqrt(sum(repmat(F(end, 2:2:end), size(combs, 1), 1).^2 .* sigma_w(2 - combs).^2, 2));
% 
% % Compute expected values (means) of y(k) for
% % each combination
% W_exp = zeros(k, 1);
% U_combined = reshape([U W_exp]', k*2, 1);
% y_mean = C*A^k*x0 + F(end, :) * U_combined;
% 
% % Calculate pdf of y
% ny = 201;
% y_range = linspace(-5, 5, ny);
% y_pdfs = normpdf(y_range, repmat(y_mean, size(combs, 1), 1), y_sigmas);
% 
% % Plot PDF
% figure(7)
% plot(y_range, sum((y_pdfs .* comb_probs), 1))
% xlabel('$y(k)$','Interpreter','Latex')
% grid on


function cd = cum_density_total(dx, dist1, dist2, p_gamma)
    cd1 = (cdf(dist1, dx) - 0.5)*2;
    cd2 = (cdf(dist2, dx) - 0.5)*2;
    cd = sum(p_gamma .* [cd1; cd2]);
end