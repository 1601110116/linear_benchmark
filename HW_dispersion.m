clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
m_min = 3;
m_max = 10;
l_num = 1;
n_num = 2;
fig_position = [50, 50, 900, 400];
line_width = 2;
marker_size = 10;
font_size = 8;
lgd_location = 'southeast';
%-----------
frqs = [23.3534, 25.8780, 26.7350, 26.5969, 26.4008, 25.8780, 25.3332, 24.8094];
gammas = [43.5309, 51.1861, 62.8595, 66.7983, 69.3115, 69.6818, 68.7938, 66.5915];

load(fullfile(code_path, 'alpha.mat'));
vte = 4.19e7 * sqrt(Tref);
nu_e = 2.91e-6 * denref * Tref ^ -1.5 * ln_lambda;
kz = 2 * pi * n_num / (height * rhos0);
c1 = kz^2 * vte^2 / (nu_e * omega_c0);
c2 = viscosity / t0 / omega_c0;
ms = m_min: m_max;
omega_r = zeros(size(ms));
omega_i = zeros(size(ms));
for m = m_min: m_max
	lambda_ml = alpha(m+1, l_num) / radius;
	b = c1 * (1+lambda_ml^2)/lambda_ml^2;
	omega_star = m / Ln^2 / (1+lambda_ml^2);
	coef1 = 1j * (b + c2*lambda_ml^2);
	coef2 = -1j * b * omega_star ...
		- lambda_ml^4 / (1+lambda_ml^2) * b * c2;
	syms omega;
	eqn = omega^2 + coef1*omega + coef2;
	res = double(solve(eqn, omega));
	if imag(res(1)) > 0
		res = res(1);
	else
		res = res(2);
	end
	res = res * omega_c0;
	omega_r(m-m_min+1) = real(res);
	omega_i(m-m_min+1) = imag(res);
end

fig = figure();
set(fig, 'Position', fig_position);

subplot(1,2,1);
plot(ms, 1e-3 * omega_r, 'b-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);
hold on;
plot(ms, 2*pi*frqs, 'r-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);
legend({'theory', 'simulation'}, 'Location', lgd_location)
xlabel('m');
ylabel('$$\omega_r\ \mathrm{\left(kHz\right)}$$', ...
	'interpreter', 'latex');
grid on;
set(gca, 'fontSize', font_size);

subplot(1,2,2);
plot(ms, 1e-3 * omega_i, 'b-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);
hold on;
plot(ms, gammas, 'r-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);
legend({'theory', 'simulation'}, 'Location', lgd_location);
xlabel('m');
ylabel('$$\gamma\ \mathrm{\left(kHz\right)}$$', ...
	'interpreter', 'latex');
grid on;
set(gca, 'fontSize', font_size);

fig_name = fullfile(fig_path, ['HW_benchmark_n', num2str(n_num), ...
	'_l', num2str(l_num), '_m', sprintf('%2.2d', m_min), ...
	'-', sprintf('%2.2d', m_max)]);
print(fig, '-dpng', fig_name);
