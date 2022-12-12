% Compute system statistics

fprintf("%15s (%1sx%1s) %7s %7s %7s\n", "Name", "w", "y", "Ts", "SNR", "N_set")


%% System 1 - SISO linear

clear

addpath("~/process-observers/")
sys_rodin_step

S = stepinfo(Gd,'SettlingTimeThreshold',0.05);
t_set = S.SettlingTime;
N_set = t_set / Ts;

Kp = dcgain(Gd);
p0 = 1;
SNR = Kp * p0 / sigma_M;

fprintf("%15s (%1dx%1d) %7g %7g %7g\n", "SISO linear", nw, ny, Ts, SNR, N_set)


%% System 2 - 2x2 linear

clear

addpath("~/process-observers/")
sys_rodin_step_2x2sym2

S = stepinfo(Gd,'SettlingTimeThreshold',0.05);
t_set = S.SettlingTime;
N_set = t_set / Ts;

Kp = dcgain(Gd);
p0 = [1; 0];
SNR = Kp(1,1) * p0(1) / sigma_M(1);

fprintf("%15s (%1dx%1d) %7g %7g %7g\n", "SISO linear", nw, ny, Ts, SNR, N_set)


%% System 3 - Grinding sim

clear

addpath("~/thesis-sims/")
rod_obs_P2DcTd4

S = stepinfo(Gd,'SettlingTimeThreshold',0.05);
t_set = S.SettlingTime;

% Or use data from grinding model

% Example 1 - step up
sim_results = readtable("~/thesis-sims/results/rod_obs_sim_1_1.csv");
figure(1); clf
make_ioplot(sim_results.Y, sim_results.t, sim_results.Pd)

dY = min(sim_results.Y) - max(sim_results.Y);
dU = max(sim_results.Pd) - min(sim_results.Pd);

nz = find(sim_results.Pd ~= 0);
k_step = nz(1);
t_step = sim_results.t(k_step);
y_95 = 0.95 * dY;
le = find(sim_results.Y <= y_95);
t_95 = sim_results.t(le(1));
t_set1 = (t_95 - t_step);

% Example 2 - step down
sim_results = readtable("~/thesis-sims/results/rod_obs_sim_1_3.csv");
figure(2); clf
make_ioplot(sim_results.Y, sim_results.t, sim_results.Pd)

dY = max(sim_results.Y) - min(sim_results.Y);
dU = min(sim_results.Pd) - max(sim_results.Pd);

t_step = 8;
y_95 = 0.95 * dY;
le = find(sim_results{sim_results.t >= t_step, 'Y'} >= -dY + y_95);
t_95 = sim_results.t(le(1)) + t_step;
t_set2 = (t_95 - t_step);

N_set = t_set2 / Ts;

SNR = abs(dY) / sigma_M(1);

fprintf("%15s (%1dx%1d) %7g %7.1f %7g\n", "SISO linear", nw, ny, Ts, SNR, N_set)