
clc
close all
clearvars

run('graphics_options.m');

%% LOAD DATA

% Use these to run the identification with different datasets
data = load("sweep_data");         % Ideal noise free data
% data = load("noisy_sweep_data"); % Noisy data (but the model is very accurate)
% data = load("wrongmodel_sweep_data");    % Noise free data (but the model NOT accurate)

%% PLOT DATA

data_plotoptions = {'LineStyle', '-', 'Color', [0      0.4471 0.7412], 'LineWidth', 2.5, 'DisplayName', 'meas'};
greybox_plotoptions  = {'LineStyle', '-', 'Color', [0.9294    0.6941    0.1255], 'LineWidth', 2, 'DisplayName', 'grey-box'};
blackbox_plotoptions  = {'LineStyle', ':', 'Color', [0.7466    0.1371    0.2981], 'LineWidth', 2, 'DisplayName', 'black-box'};

figure;
tiledlayout(2, 1, 'TileSpacing', 'tight');

tiles(1) = nexttile; hold on; grid on;
ylabel('$\ddot{z}$ [mm/s$^2$]');
plot(data.time, data.y5*1000, data_plotoptions{:});

tiles(2) = nexttile; hold on; grid on;
ylabel('$F$ [N]');
xlabel('$time$ [s]');
plot(data.time, data.u*1000, data_plotoptions{:});

linkaxes(tiles, 'x');
tiles(1).XTickLabel = '';
xlim([data.time(1) data.time(end)]);
drawnow;

%% ANALYZE DATA IN FREQUENCY DOMAIN

f_patch_range = [0.9 15];
f_plot_range = [10^-1 10^2];

% function compute_io_spectrum is used to compute the raw spectrum using
% fft and a filtered spectrum using a filtered version of fft based on
% pwelch power spectral density computation
[raw_spectrum_1, filtered_spectrum_1] = compute_io_spectrum(data.y1, data.u, 1/data.Ts);
[raw_spectrum_2, filtered_spectrum_2] = compute_io_spectrum(data.y2, data.u, 1/data.Ts);
[raw_spectrum_3, filtered_spectrum_3] = compute_io_spectrum(data.y3, data.u, 1/data.Ts);
[raw_spectrum_4, filtered_spectrum_4] = compute_io_spectrum(data.y4, data.u, 1/data.Ts);
[raw_spectrum_5, filtered_spectrum_5] = compute_io_spectrum(data.y5, data.u, 1/data.Ts);

filt_spectrum_plotoptions = {'LineStyle', '-', 'Color', [0      0.4471 0.7412], 'LineWidth', 2.5, 'DisplayName', 'filt spectrum'};
raw_spectrum_plotoptions  = {'LineStyle', '-', 'Color', [0.5000 0.7235 0.8706 0.5], 'LineWidth', 2.5, 'DisplayName', 'raw spectrum'};
patch_plotoptions = {[1 1 1]*0.85, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off'};

clearvars tiles
figure;
tiledlayout(2, 1, 'TileSpacing', 'compact');

tiles(1) = nexttile; grid on; hold on;
ylabel('$|G(j \omega)|$ [dB]');
plot(raw_spectrum_1.freq, 20*log10(raw_spectrum_1.mag), raw_spectrum_plotoptions{:}); hold on; 
plot(filtered_spectrum_1.freq, 20*log10(filtered_spectrum_1.mag), filt_spectrum_plotoptions{:});

tiles(2) = nexttile; hold on; grid on;
ylabel('$\angle{G(j \omega)}$ [deg]');
xlabel('$freq$ [Hz]');
plot(raw_spectrum_1.freq, raw_spectrum_1.phase*180/pi, raw_spectrum_plotoptions{:});
plot(filtered_spectrum_1.freq, filtered_spectrum_1.phase*180/pi, filt_spectrum_plotoptions{:});

mag_plot_range = [min(20*log10(raw_spectrum_1.mag)) max(20*log10(raw_spectrum_1.mag))];
ph_plot_range = [min(raw_spectrum_1.phase(filtered_spectrum_1.freq<f_patch_range(2))*180/pi) max(raw_spectrum_1.phase(filtered_spectrum_1.freq<f_patch_range(2))*180/pi)];
mag_plot_range = mag_plot_range + 0.1*[-1 +1]*diff(mag_plot_range);
ph_plot_range = ph_plot_range + 0.1*[-1 +1]*diff(ph_plot_range);

linkaxes(tiles, 'x');
tiles(1).XTickLabel = '';
xlim(f_plot_range);
tiles(1).XScale = 'log';
tiles(2).XScale = 'log';
tiles(1).YLim = mag_plot_range;
tiles(2).YLim = ph_plot_range;
legend(tiles(1));

patch(tiles(1), [f_plot_range(1); f_patch_range(1); f_patch_range(1); f_plot_range(1)], [mag_plot_range(1); mag_plot_range(1); mag_plot_range(2); mag_plot_range(2)], patch_plotoptions{:});
patch(tiles(1), [f_patch_range(2); f_plot_range(2); f_plot_range(2); f_patch_range(2)], [mag_plot_range(1); mag_plot_range(1); mag_plot_range(2); mag_plot_range(2)], patch_plotoptions{:});
patch(tiles(2), [f_plot_range(1); f_patch_range(1); f_patch_range(1); f_plot_range(1)], [ph_plot_range(1); ph_plot_range(1); ph_plot_range(2); ph_plot_range(2)], patch_plotoptions{:});
patch(tiles(2), [f_patch_range(2); f_plot_range(2); f_plot_range(2); f_patch_range(2)], [ph_plot_range(1); ph_plot_range(1); ph_plot_range(2); ph_plot_range(2)], patch_plotoptions{:});

title(tiles(1), 'TF from $F$ to $\dot{z}$');
drawnow;

%% GREY-BOX IDENTIFICATION

% Known parameters:
mt = 64.069558;
kt = 200000;

% theta = [M, k, c]
% Compute an initial guess theta_0 and upper and lower bounds for the
% optimization
theta_ub = [900 250000 15000]';
theta_lb = [400 30000 2000]';
theta_0 = [600 10000 5000]';

% Cost function definition
my_cost_fcn_handle = @(theta) my_cost_fcn(theta, data, mt, kt);

% Run optimization
% theta_opt = particleswarm(my_cost_fcn_handle, 3, theta_lb, theta_ub, optimoptions('particleswarm', 'Display', 'iter'));
theta_opt = fmincon(my_cost_fcn_handle, theta_0, [], [], [], [], theta_lb, theta_ub, [], optimoptions('fmincon', 'Display', 'iter'));

% Get optimal values
M_id = theta_opt(1);
k_id = theta_opt(2);
c_id = theta_opt(3);

% Compute grey-box system model
greybox_sys = build_QC_model(M_id, k_id, c_id, mt, kt);
greybox_tf = tf(greybox_sys);
greybox_y = lsim(greybox_sys, data.u, data.time);

%% BLACK-BOX (time domain) IDENTIFICATION
% We will use tfest function, which is a very powerful function able to
% estimate a model of a given order from data both in time domain and
% frequency domain (in this case we are using the time domain version)

% Choose the number of poles of each transfer function
n_poles = [4 4 4 4 4]';
n_zeros = [3 2 3 2 4]';
delay   = [0 0 0 0 0]';

% Identify the 5 tfs using tfest
id_data = iddata([data.y1, data.y2, data.y3, data.y4, data.y5], data.u, data.Ts);
blackbox_sys = tfest(id_data, n_poles, n_zeros, delay);
blackbox_y = lsim(blackbox_sys, data.u, data.time);

%% TF COMPARISON

fprintf(1, "\n\nGrey-box TF from u to y2:\n");
zpk(minreal(greybox_tf(2)))
fprintf(1, "\n\nGrey-box TF from u to y5:\n");
zpk(minreal(greybox_tf(5)))
fprintf(1, "Black-box TF from u to y2:\n");
zpk(minreal(tf(blackbox_sys(2))))
fprintf(1, "Black-box TF from u to y5:\n");
zpk(minreal(tf(blackbox_sys(5))))

%% PLOTS in time domain

clearvars tiles;
figure;
tiledlayout(5, 1, 'TileSpacing', 'tight');

tiles(1) = nexttile; hold on; grid on;
ylabel('$\dot{z}$ [mm/s]');
plot(data.time, data.y1*1000, data_plotoptions{:});
plot(data.time, greybox_y(:, 1)*1000, greybox_plotoptions{:});
plot(data.time, blackbox_y(:, 1)*1000, blackbox_plotoptions{:});

tiles(2) = nexttile; hold on; grid on;
ylabel('${z}$ [mm]');
plot(data.time, data.y2*1000, data_plotoptions{:});
plot(data.time, greybox_y(:, 2)*1000, greybox_plotoptions{:});
plot(data.time, blackbox_y(:, 2)*1000, blackbox_plotoptions{:});

tiles(3) = nexttile; hold on; grid on;
ylabel('$\dot{z_t}$ [mm/s]');
plot(data.time, data.y3*1000, data_plotoptions{:});
plot(data.time, greybox_y(:, 3)*1000, greybox_plotoptions{:});
plot(data.time, blackbox_y(:, 3)*1000, blackbox_plotoptions{:});

tiles(4) = nexttile; hold on; grid on;
ylabel('${z_t}$ [mm]');
plot(data.time, data.y4*1000, data_plotoptions{:});
plot(data.time, greybox_y(:, 4)*1000, greybox_plotoptions{:});
plot(data.time, blackbox_y(:, 4)*1000, blackbox_plotoptions{:});

tiles(5) = nexttile; hold on; grid on;
ylabel('$\ddot{z}$ [mm/s$^2$]');
xlabel('$time$ [s]');
plot(data.time, data.y5*1000, data_plotoptions{:});
plot(data.time, greybox_y(:, 5)*1000, greybox_plotoptions{:});
plot(data.time, blackbox_y(:, 5)*1000, blackbox_plotoptions{:});

linkaxes(tiles, 'x');
tiles(1).XTickLabel = '';
tiles(2).XTickLabel = '';
tiles(3).XTickLabel = '';
tiles(4).XTickLabel = '';

legend(tiles(1), 'Location', 'northeastoutside');

xlim([data.time(1) data.time(end)]);
drawnow;

%% FREQUENCY DOMAIN ANALYSIS

omega_vec = filtered_spectrum_1.freq*2*pi;

% Grey-Box identified model:
[ gb_mag1, gb_phase1 ] = bode(greybox_tf(1,1).num, greybox_tf(1,1).den, omega_vec);
[ gb_mag2, gb_phase2 ] = bode(greybox_tf(2,1).num, greybox_tf(2,1).den, omega_vec);
[ gb_mag3, gb_phase3 ] = bode(greybox_tf(3,1).num, greybox_tf(3,1).den, omega_vec);
[ gb_mag4, gb_phase4 ] = bode(greybox_tf(4,1).num, greybox_tf(4,1).den, omega_vec);
[ gb_mag5, gb_phase5 ] = bode(greybox_tf(5,1).num, greybox_tf(5,1).den, omega_vec);

gb_phase1 = wrapToPi(gb_phase1*pi/180 - filtered_spectrum_1.phase)*180/pi + filtered_spectrum_1.phase*180/pi;
gb_phase2 = wrapToPi(gb_phase2*pi/180 - filtered_spectrum_2.phase)*180/pi + filtered_spectrum_2.phase*180/pi;
gb_phase3 = wrapToPi(gb_phase3*pi/180 - filtered_spectrum_3.phase)*180/pi + filtered_spectrum_3.phase*180/pi;
gb_phase4 = wrapToPi(gb_phase4*pi/180 - filtered_spectrum_4.phase)*180/pi + filtered_spectrum_4.phase*180/pi;
gb_phase5 = wrapToPi(gb_phase5*pi/180 - filtered_spectrum_5.phase)*180/pi + filtered_spectrum_5.phase*180/pi;

% Black-Box identified model:
[ bb_mag1, bb_phase1 ] = bode(blackbox_sys(1,1).num, blackbox_sys(1,1).den, omega_vec);
[ bb_mag2, bb_phase2 ] = bode(blackbox_sys(2,1).num, blackbox_sys(2,1).den, omega_vec);
[ bb_mag3, bb_phase3 ] = bode(blackbox_sys(3,1).num, blackbox_sys(3,1).den, omega_vec);
[ bb_mag4, bb_phase4 ] = bode(blackbox_sys(4,1).num, blackbox_sys(4,1).den, omega_vec);
[ bb_mag5, bb_phase5 ] = bode(blackbox_sys(5,1).num, blackbox_sys(5,1).den, omega_vec);

bb_phase1 = wrapToPi(bb_phase1*pi/180 - filtered_spectrum_1.phase)*180/pi + filtered_spectrum_1.phase*180/pi;
bb_phase2 = wrapToPi(bb_phase2*pi/180 - filtered_spectrum_2.phase)*180/pi + filtered_spectrum_2.phase*180/pi;
bb_phase3 = wrapToPi(bb_phase3*pi/180 - filtered_spectrum_3.phase)*180/pi + filtered_spectrum_3.phase*180/pi;
bb_phase4 = wrapToPi(bb_phase4*pi/180 - filtered_spectrum_4.phase)*180/pi + filtered_spectrum_4.phase*180/pi;
bb_phase5 = wrapToPi(bb_phase5*pi/180 - filtered_spectrum_5.phase)*180/pi + filtered_spectrum_5.phase*180/pi;

clearvars tiles
figure;
tiledlayout(2, 1, 'TileSpacing', 'compact');

tiles(1) = nexttile; grid on; hold on;
ylabel('$|G(j \omega)|$ [dB]');
plot(raw_spectrum_1.freq, 20*log10(raw_spectrum_1.mag), raw_spectrum_plotoptions{:}); hold on; 
plot(filtered_spectrum_1.freq, 20*log10(filtered_spectrum_1.mag), filt_spectrum_plotoptions{:});
plot(omega_vec/2/pi, 20*log10(gb_mag1), greybox_plotoptions{:});
plot(omega_vec/2/pi, 20*log10(bb_mag1), blackbox_plotoptions{:});

tiles(2) = nexttile; hold on; grid on;
ylabel('$\angle{G(j \omega)}$ [deg]');
xlabel('$freq$ [Hz]');
plot(raw_spectrum_1.freq, raw_spectrum_1.phase*180/pi, raw_spectrum_plotoptions{:});
plot(filtered_spectrum_1.freq, filtered_spectrum_1.phase*180/pi, filt_spectrum_plotoptions{:});
plot(omega_vec/2/pi, gb_phase1, greybox_plotoptions{:});
plot(omega_vec/2/pi, bb_phase1, blackbox_plotoptions{:});

mag_plot_range = [min(20*log10(raw_spectrum_1.mag)) max(20*log10(raw_spectrum_1.mag))];
ph_plot_range = [min(raw_spectrum_1.phase(filtered_spectrum_1.freq<f_patch_range(2))*180/pi) max(raw_spectrum_1.phase(filtered_spectrum_1.freq<f_patch_range(2))*180/pi)];
mag_plot_range = mag_plot_range + 0.1*[-1 +1]*diff(mag_plot_range);
ph_plot_range = ph_plot_range + 0.1*[-1 +1]*diff(ph_plot_range);

linkaxes(tiles, 'x');
tiles(1).XTickLabel = '';
xlim(f_plot_range);
tiles(1).XScale = 'log';
tiles(2).XScale = 'log';
tiles(1).YLim = mag_plot_range;
tiles(2).YLim = ph_plot_range;
legend(tiles(1));

patch(tiles(1), [f_plot_range(1); f_patch_range(1); f_patch_range(1); f_plot_range(1)], [mag_plot_range(1); mag_plot_range(1); mag_plot_range(2); mag_plot_range(2)], patch_plotoptions{:});
patch(tiles(1), [f_patch_range(2); f_plot_range(2); f_plot_range(2); f_patch_range(2)], [mag_plot_range(1); mag_plot_range(1); mag_plot_range(2); mag_plot_range(2)], patch_plotoptions{:});
patch(tiles(2), [f_plot_range(1); f_patch_range(1); f_patch_range(1); f_plot_range(1)], [ph_plot_range(1); ph_plot_range(1); ph_plot_range(2); ph_plot_range(2)], patch_plotoptions{:});
patch(tiles(2), [f_patch_range(2); f_plot_range(2); f_plot_range(2); f_patch_range(2)], [ph_plot_range(1); ph_plot_range(1); ph_plot_range(2); ph_plot_range(2)], patch_plotoptions{:});

title(tiles(1), 'TF from $F$ to $\dot{z}$');
drawnow;

%% FUNCTIONS

function J_cost = my_cost_fcn( theta, data, mt, kt )    
    
    % Theta is the vector of parameters:
    M = theta(1);
    k = theta(2);
    c = theta(3);
    
    % Get the quarter car model in discrete time
    sys = build_QC_model(M, k, c, mt, kt, data.Ts);
    
    % Simulate the QC model with the current M, k, c
    sim_y = lsim(sys, data.u, data.time);
    
    % Compute the root mean square difference between each simulated output
    % and the real data. We divide each term by the data standard deviation
    % as normalization term.
    J_cost_1 = rms(sim_y(:,1)-data.y1)/std(data.y1);
    J_cost_2 = rms(sim_y(:,2)-data.y2)/std(data.y2);
    J_cost_3 = rms(sim_y(:,3)-data.y3)/std(data.y3);
    J_cost_4 = rms(sim_y(:,4)-data.y4)/std(data.y4);
    J_cost_5 = rms(sim_y(:,5)-data.y5)/std(data.y5); 
    
    % Compute the overall cost of the current experiment
    J_cost = J_cost_1 + J_cost_2 + J_cost_3 + J_cost_4 + J_cost_5;
                
end

% Function to build the Quarter-Car (QC) model in continuous or discrete
% time, depending on the presence of the last parameter Ts
function sys = build_QC_model(M, k, c, mt, kt, Ts)
    
    % Quarter-car model:  u = F
    %                     y = [zdot z ztdot zt zdotdot]
    %                     x = [zdot z ztdot zt]
    
    % System Matrices in continuous time
    A = [-c/M   -k/M     c/M          k/M;
            1      0       0            0;
         c/mt   k/mt   -c/mt   -(k+kt)/mt;
            0      0       1            0];
    B = [  1/M;
             0;
         -1/mt;
             0];
    C = [eye(4);
         -c/M   -k/M   c/M   k/M];
    D = [zeros(4,1); 1/M];
    
    
    if exist('Ts', 'var')
        % Create continuous-time system
        sys_c = ss(A, B, C, D);

        % Discretize CT system in according to sampling time Ts
        sys = c2d(sys_c, Ts, 'tustin');
    else
        % Create continuous-time system
        sys = ss(A, B, C, D);
    end
    
end

function [raw_spectrum, filtered_spectrum] = compute_io_spectrum(y_data, u_data, f_sampling)
    
%% RAW SPECTRUM
    
    % We will compute the Fourier transform of the input and output signals
    % and compute the spectrum as the simple ratio between output and input
    % spectra    
    N = length(y_data);    
    freq = (f_sampling*(0:floor(N/2))/N)';
    fftout = fft(y_data, N);
    fftin  = fft(u_data, N);
    
    modY                 = nan(size(freq));
    modY(1)              = abs(fftout(1))/N; % The first element is normalized by N
    modY(2:floor(N/2)+1) = abs(fftout(2:floor(N/2)+1))*2/N;    
    phaseY               = angle(fftout(1:floor(N/2)+1));
    
    modU                 = nan(size(freq));
    modU(1)              = abs(fftin(1))/N; % The first element is normalized by N
    modU(2:floor(N/2)+1) = abs(fftin(2:floor(N/2)+1))*2/N;    
    phaseU               = angle(fftin(1:floor(N/2)+1));
    
    % Compute raw tf spectrum as ratio between output and input spectra    
    raw_spectrum.freq  = freq;
    raw_spectrum.mag   = modY./modU;
    raw_spectrum.phase = unwrap(wrapToPi(phaseY - phaseU));
    
%% FILTERED SPECTRUM
    
    % To compute the filtered spectrum we can use the function tfestimate,
    % which uses a short-time fourier transform of the signals. It computes
    % the spectrum as out_spectruum/in_spectrum but the out and in spectra
    % are computed using Welch's averaged periodogram method.
    spectrum = tfestimate(u_data, y_data, [], [], freq, f_sampling);
    filtered_spectrum.freq  = freq;
    filtered_spectrum.mag   = reshape(abs(spectrum), [], 1);
    filtered_spectrum.phase = reshape(unwrap(wrapToPi(angle(spectrum))), [], 1);
    
    raw_spectrum.phase   = wrapToPi(raw_spectrum.phase - filtered_spectrum.phase) + filtered_spectrum.phase;
    
end

