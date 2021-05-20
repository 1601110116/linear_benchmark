clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
l_num = 1;
fig_position = [50, 50, 900, 400];
font_size = 8;
lgd_location = 'southeast';
nu_e = 2.91e-6 * denref * Tref ^ -1.5 * ln_lambda;
c2 = 8 * viscosity / t0 / omega_c0;
Ln = 3 / rhos0;
m_max = 12;
n_max = 3;  % do not change
line_width = 2;
marker_size = 10;
%-----------

load(fullfile(code_path, 'alpha.mat'));
vte = 4.19e7 * sqrt(Tref);
omega_r = zeros(m_max, n_max);
omega_i = zeros(m_max, n_max);
for n = 1: n_max
	kz = 2 * pi * n / (height * rhos0);
	c1 = kz^2 * vte^2 / (nu_e * omega_c0);
	for m = 1: m_max
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
		omega_r(m, n) = real(res);
		omega_i(m, n) = imag(res);
	end
end

fig = figure();
set(fig, 'Position', fig_position);

ms = 1: m_max;
subplot(1,2,1);
plot(ms, 1e-3 * omega_r(:, 1), 'b-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);  hold on;
plot(ms, 1e-3 * omega_r(:, 2), 'r-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);  
plot(ms, 1e-3 * omega_r(:, 3), 'k-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);  
legend({'n=1', 'n=2', 'n=3'}, 'Location', lgd_location);
xlabel('m');
ylabel('$$\omega_r\ \mathrm{\left(kHz\right)}$$', ...
	'interpreter', 'latex');
grid on;
set(gca, 'fontSize', font_size);

subplot(1,2,2);
plot(ms, 1e-3 * omega_i(:, 1), 'b-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size);  hold on;
plot(ms, 1e-3 * omega_i(:, 2), 'r-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size); 
plot(ms, 1e-3 * omega_i(:, 3), 'k-o', 'lineWidth', line_width, ...
	'MarkerSize', marker_size); 
legend({'n=1', 'n=2', 'n=3'}, 'Location', lgd_location);
xlabel('m');
ylabel('$$\gamma\ \mathrm{\left(kHz\right)}$$', ...
	'interpreter', 'latex');
grid on;
set(gca, 'fontSize', font_size);

fig_name = fullfile(fig_path, 'HW_theory');
print(fig, '-dpng', fig_name);
