clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
start_diag = 30;
end_diag = 70; %last_diagose;
m_num = 8;
n_num = 2;
l_max = 30;
zdiag = length(z) / 2;
ydiag = length(x) / 2;
fig_position = [50, 50, 800, 600];
font_size = 15;
%-----------

build_grid_2d;
HW_mn_data = fullfile('data', ...
	['HW_m', sprintf('%2.2d', m_num), 'n', num2str(n_num)]);
HW_mn_diag = fullfile('diagnose', ...
	['HW_m', sprintf('%2.2d', m_num), 'n', num2str(n_num)]);
if ~exist(HW_mn_data, 'dir')
	mkdir(HW_mn_data);
end
if ~exist(HW_mn_diag, 'dir')
	mkdir(HW_mn_diag);
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
	den_m = filter_m(den_n, m_num, l_max);
	den_m = den_m ./ den0;
	phi_n = filter_z(phi, n_num);
	% Te=Tref is assumed
	phi_m = filter_m(phi_n, m_num, l_max);
	fig = figure('Visible', 'off');
	set(fig, 'Position', fig_position);

	subplot(2,2,1);
	pcolor(rhos0*y2d, rhos0*x2d, den_m(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	xlabel('y (cm)');
	ylabel('x (cm)');
	title('HW $$\tilde{n}$$', 'interpreter', 'latex');
	set(gca, 'fontSize', font_size)

	subplot(2,2,2);
	pcolor(rhos0*y2d, rhos0*x2d, phi_m(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	xlabel('y (cm)');
	ylabel('x (cm)');
	title('HW $$\tilde{\phi}$$', 'interpreter', 'latex');
	set(gca, 'fontSize', font_size)

	subplot(2,2,3);
	pcolor(zX, xX, squeeze(den_m(:, ydiag, :)));
	colormap jet;  colorbar;  shading interp;
	xlabel('z (cm)');
	ylabel('x (cm)');
	set(gca, 'fontSize', font_size)

	subplot(2,2,4);
	pcolor(zX, xX, squeeze(phi_m(:, ydiag, :)));
	colormap jet;  colorbar;  shading interp;
	xlabel('z (cm)');
	ylabel('x (cm)');
	set(gca, 'fontSize', font_size)

	file_name = ['m', sprintf('%2.2d', m_num), 'n', num2str(n_num), ...
		'z', sprintf('%2.2d', zdiag), 't', sprintf('%4.4d', idiag)];
	save(fullfile(HW_mn_data, file_name), 'den_m', 'phi_m');
	print(fig, '-dpng', fullfile(HW_mn_diag, file_name));
end


