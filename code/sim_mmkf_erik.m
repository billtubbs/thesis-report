% Simulate AFMM procedure

addpath('~/process-observers/')

% Set random number seed
rng(0)

% Simulation parameters
Ts = 1; nT = 15;
t = Ts*(0:nT)';

% RODD random variable parameters
epsilon = 0.01;
sigma_wp = [0.01 1];

% Generate random shock signal
Gamma = zeros(nT+1, 1);
Gamma(t == 4) = 1;
Wp = 0 .* randn(nT+1, 1);
Wp(Gamma == 1) = 1 .* sigma_wp(2);

% Step Disturbance Process
Hd = rodd_tf(1, 1, 1, 1);
P = lsim(Hd,Wp,t);

figure(1); clf
stairs(t, P, 'Linewidth', 2)
grid on
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$p(k)$', 'Interpreter', 'latex')

% Process noise standard deviation
sigma_W = [0; 0];

% Discrete time state space model
A = [0.7 1;
     0 1];
B = [1 0;
     0 1];
C = [0.3 0];
D = zeros(1, 2);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Prepare inputs to simulation
U = zeros(size(t));
%U(t >= 1) = -1;
U_sim = [U Wp];

% Default initial condition
x0 = zeros(n, 1);
p0 = 0;

% Check simulation output is correct
[Y, t, X] = lsim(Gpss,U_sim,t);

% Add measurement noise for plant
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));

% Designate which input and output variables are
% measured
u_meas = [true; false];
y_meas = true;

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Measurement noise standard deviation
sigma_M = 0.001;

% Multiple model observer with AFMM algorithm
P0 = 0.01*eye(n);
Q0 = diag([0.1 0]);
R = sigma_M^2;
f = nT+1;  % sequence history length
n_filt = 5;  % number of filters
n_min = 2;  % minimum life of cloned filters
obs = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,'MKF_AFMM');

% Prepare simulation varibles
k = (0:nT)';
X_est = nan(nT+1,n);
Y_est = nan(nT+1,ny);
E_obs = nan(nT+1,ny);

% Arrays to store observer variables
switch obs.label
    case {'MKF'}
        n_filt = obs.n_filt;
        MKF_p_seq_g_Yk = nan(nT+1, n_filt);
    case {'MKF_AFMM'}
        n_filt = obs.n_filt;
        MKF_p_seq_g_Yk = nan(nT+1, n_filt);
        AFMM_f_main = nan(nT+1, numel(obs.f_main));
        AFMM_f_hold = nan(nT+1, numel(obs.f_hold));
    otherwise
        n_filt = 1;
end
K_obs = cell(nT+1, n_filt);
trP_obs = cell(nT+1, n_filt);

% Start simulation at k = 0
for i = 1:nT+1

    % For debugging:
    %fprintf("t = %f\n", t(i));

    % Process measurements
    uk_m = U(i,:)';
    yk_m = Y_m(i,:)';

    % Record observer estimates and output errors
    X_est(i, :) = obs.xkp1_est';
    Y_est(i, :) = obs.ykp1_est';
    E_obs(i, :) = yk_m' - obs.ykp1_est';

    % Kalman update equations
    % Update observer gains and covariance matrix
    switch obs.type

        case {'MKF', 'MKF_RODD'}
            obs = update_MKF(obs, uk_m, yk_m);

            % Record filter gains and covariance matrices
            for j=1:obs.n_filt
                K_obs{i, j} = obs.filters{j}.K';
                trP_obs{i, j} = trace(obs.filters{j}.P);
            end

            % Record filter conditional probabilities
            MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

        case 'MKF_AFMM'

            obs = update_AFMM(obs, uk_m, yk_m);

            % Record filter arrangement before update
            AFMM_f_main(i, :) = obs.f_main;
            AFMM_f_hold(i, :) = obs.f_hold;

            % Record filter conditional probabilities
            MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

            % Record filter gains and covariance matrices
            for j=1:obs.n_filt
                K_obs{i, j} = obs.filters{j}.K';
                trP_obs{i, j} = trace(obs.filters{j}.P);
            end

        otherwise
            error('Observer type not valid')
        disp('pause')
    end

end

% Display results table
table(k, t, Gamma, P, Wp, U, X, Y)

figure(2); clf

ax1 = subplot(3,1,1);
plot(t, Y_m, 'o'); hold on
plot(t, Y, '-', 'Linewidth', 2);
plot(t, Y_est, '-', 'Linewidth', 2);
grid on
legend('$y_m(k)$', '$y(k)$', '$\hat{y}(k)$', 'Interpreter', 'latex', 'location' , 'best')

ax2 = subplot(3,1,2);
plot(t, X_est, 'Linewidth', 2); hold on
grid on
legend('$\hat{x}_1(k)$', '$\hat{x}_2(k)$', 'Interpreter', 'latex', 'location' , 'best')

ax3 = subplot(3,1,3);
stairs(t, [U P], 'Linewidth', 2); hold on
grid on
legend('$u(k)$', '$p(k)$', 'Interpreter', 'latex', 'location' , 'best')
xlabel('$t$', 'Interpreter', 'latex')

linkaxes([ax1 ax2 ax3], 'x')

% Display AFMM filter groupings and probabilities
switch obs.label
    case {'MKF_AFMM'}
    f_hold = AFMM_f_hold;
    f_main = AFMM_f_main;
    p_seq_g_Yk = round(MKF_p_seq_g_Yk, 3);
    [table(k, t) array2table(f_hold) array2table(f_main) array2table(p_seq_g_Yk)]
end