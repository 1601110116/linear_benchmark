clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
start_diag = 30;
end_diag = 31; %last_diagose;
m_num = 8;
n_num = 2;
l_num = 1;
zdiag = length(z) / 2;
ydiag = length(x) / 2;
fig_position = [50, 50, 800, 600];
font_size = 15;
x_ind = 30;
y_ind = ydiag;
marker_size = 10;
%-----------

build_grid_2d;
HW_mnl_data = fullfile('data', ...
	['HW_m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
	'l', num2str(l_num)]);
HW_mnl_diag = fullfile('diagnose', ...
	['HW_m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
	'l', num2str(l_num)]);
if ~exist(HW_mnl_data, 'dir')
	mkdir(HW_mnl_data);
end
if ~exist(HW_mnl_diag, 'dir')
	mkdir(HW_mnl_diag);
end

% calculate the n_0 for the new normalization
den0 = N0 * exp(-0.5 * (r2d / Ln) .^ 2);

for idiag = start_diag: end_diag
	close;
	disp(['visualizing HW den and phi, step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
    load(data_file, 'den', 'phi');

	den_n = filter_z(den, n_num);
	den_n = den_n ./ den0;
	[~, den_ml] = fbt_ml(den_n, m_num, l_num);

	phi_n = filter_z(phi, n_num);
	% Te=Tref is assumed
	[~, phi_ml] = fbt_ml(phi_n, m_num, l_num);

	fig = figure('Visible', 'off');
	set(fig, 'Position', fig_position);

	subplot(2,2,1);
	pcolor(rhos0*y2d, rhos0*x2d, den_ml(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	hold on;
	plot(rhos0*y2d(x_ind, y_ind), rhos0*x2d(x_ind, y_ind), 'w+', ...
		'MarkerSize', marker_size);  hold off;
	xlabel('y (cm)');
	ylabel('x (cm)');
	title('HW $$\tilde{n}$$', 'interpreter', 'latex');
	set(gca, 'fontSize', font_size)

	subplot(2,2,2);
	pcolor(rhos0*y2d, rhos0*x2d, phi_ml(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	hold on;
	plot(rhos0*y2d(x_ind, y_ind), rhos0*x2d(x_ind, y_ind), 'w+', ...
		'MarkerSize', marker_size);  hold off;
	xlabel('y (cm)');
	ylabel('x (cm)');
	title('HW $$\tilde{\phi}$$', 'interpreter', 'latex');
	set(gca, 'fontSize', font_size)

	subplot(2,2,3);
	pcolor(zX, xX, squeeze(den_ml(:, ydiag, :)));
	colormap jet;  colorbar;  shading interp;
	hold on;
	plot(zX(zdiag), xX(x_ind), 'w+', 'MarkerSize', marker_size);  hold off;
	xlabel('z (cm)');
	ylabel('x (cm)');
	set(gca, 'fontSize', font_size)

	subplot(2,2,4);
	pcolor(zX, xX, squeeze(phi_ml(:, ydiag, :)));
	colormap jet;  colorbar;  shading interp;
	hold on;
	plot(zX(zdiag), xX(x_ind), 'w+', 'MarkerSize', marker_size);  hold off;
	xlabel('z (cm)');
	ylabel('x (cm)');
	set(gca, 'fontSize', font_size)

	file_name = ['m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
		'l', num2str(l_num), ...
		'z', sprintf('%2.2d', zdiag), 't', sprintf('%4.4d', idiag)];
	save(fullfile(HW_mnl_data, file_name), 'den_ml', 'phi_ml');
	print(fig, '-dpng', fullfile(HW_mnl_diag, file_name));
end


