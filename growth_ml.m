clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
start_diag = 250;
end_diag = 800;
m_num = 10;
n_num = 2;
l_num = 1;
zdiag = length(z) / 2;
ydiag = length(x) / 2;
x_ind = 30;
y_ind = ydiag;
z_ind = zdiag;
fig_position = [50, 50, 1600, 400];
font_size = 15;
line_width = 2;
marker_size = 12;
%-----------

build_grid_2d;
HW_mnl_data = fullfile('data', ...
    ['HW_m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
    'l', num2str(l_num)]);
HW_mnl_diag = fullfile('diagnose', ...
    ['HW_m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
	'l', num2str(l_num)]);

ndiag = end_diag - start_diag + 1;
den_t = zeros(1, ndiag);
phi_t = zeros(1, ndiag);
for idiag = start_diag: end_diag
	file_name = ['m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
        'l', num2str(l_num), ... 
        'z', sprintf('%2.2d', zdiag), 't', sprintf('%4.4d', idiag)];
	data_file = fullfile(HW_mnl_data, file_name);
	load(data_file, 'den_ml', 'phi_ml');
	den_t(idiag-start_diag+1) = den_ml(x_ind, y_ind, z_ind);
	phi_t(idiag-start_diag+1) = phi_ml(x_ind, y_ind, z_ind);
end

tX_ms = 1e3 * t0 * dt * nt_per_diagnose * (start_diag: end_diag);
ln_den_t = log(abs(den_t));
ln_phi_t = log(abs(phi_t));

den_peaks = find(diff(sign(diff(ln_den_t))) == -2) + 1;
% in ms
den_peaks_time = tX_ms(den_peaks);
den_peaks_lnden = ln_den_t(den_peaks);
den_fit_poly = polyfit(den_peaks_time, den_peaks_lnden, 1);
den_periods = den_peaks_time(3:end) - den_peaks_time(1:end-2);
den_frqs = 1 ./ den_periods;
den_mean_frq = mean(den_frqs);

phi_peaks = find(diff(sign(diff(ln_phi_t))) == -2) + 1;
% in ms
phi_peaks_time = tX_ms(phi_peaks);
phi_peaks_lnphi = ln_phi_t(phi_peaks);
phi_fit_poly = polyfit(phi_peaks_time, phi_peaks_lnphi, 1);
phi_periods = phi_peaks_time(3:end) - phi_peaks_time(1:end-2);
phi_frqs = 1 ./ phi_periods;
phi_mean_frq = mean(phi_frqs);

fig = figure();
set(fig, 'Position', fig_position);

subplot(1,2,1);
yyaxis left;
plot(tX_ms, ln_den_t, '--', 'lineWidth', line_width);  hold on;
plot(tX_ms, den_fit_poly(1)*tX_ms + den_fit_poly(2), ...
	'-', 'lineWidth', line_width);
plot(den_peaks_time, den_peaks_lnden, '+', 'MarkerSize', marker_size);
grid on;
xlabel('t (ms)');
ylabel('$$\ln\left|\tilde{n}\right|$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size);
y_lim = get(gca, 'YLim');
yyaxis right;
plot(den_peaks_time(2:end-1), den_frqs, '-', 'lineWidth', line_width);
yline(den_mean_frq, '--', 'lineWidth', line_width);
ylabel('$$f\ \mathrm{\left(kHz\right)}$$', 'interpreter', 'latex');
title(['$$f=', num2str(den_mean_frq), '$$ kHz, $$\gamma=', ...
	num2str(den_fit_poly(1)), '$$ kHz'], 'interpreter', 'latex');
set(gca, 'FontSize', font_size);

subplot(1,2,2);
yyaxis left;
plot(tX_ms, ln_phi_t, '--', 'lineWidth', line_width);  hold on;
plot(tX_ms, phi_fit_poly(1)*tX_ms + phi_fit_poly(2), ...
	'-', 'lineWidth', line_width);
plot(phi_peaks_time, phi_peaks_lnphi, '+', 'MarkerSize', marker_size);
grid on;
xlabel('t (ms)');
ylabel('$$\ln\left|\tilde{\phi}\right|$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size);
set(gca, 'YLim', y_lim);
yyaxis right;
plot(phi_peaks_time(2:end-1), phi_frqs, '-', 'lineWidth', line_width);
yline(phi_mean_frq, '--', 'lineWidth', line_width);
ylabel('$$f\ \mathrm{\left(kHz\right)}$$', 'interpreter', 'latex');
title(['$$f=', num2str(phi_mean_frq), '$$ kHz, $$\gamma=', ...
	num2str(phi_fit_poly(1)), '$$ kHz'], 'interpreter', 'latex');
set(gca, 'FontSize', font_size);

fig_name = fullfile(HW_mnl_diag, ['linear_x', num2str(x_ind), ...
	'_y', num2str(y_ind), '_z', num2str(z_ind)]);
print(fig, '-dpng', fig_name);

