clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
start_diag = 0;
end_diag = last_diagnose;
l_max = 30;
line_width = 2;
fig_position = [50, 50, 300, 300];
font_size = 15;
%-----------

build_grid;
build_grid_2d;
for idiag = start_diag: end_diag
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'den');
	denp = field2pol(den);
	den0_interp = denref * squeeze(zonal_average(denp));
	den0_fbt = fbt_mr(den, r, 0, l_max);
	den0_fbt = denref * mean(den0_fbt, 3);
	lnn0_interp = log(den0_interp);
	lnn0_fbt = log(den0_fbt);
	rX = rhos0 * r;
	fig = figure('Visible', 'off');
	plot(rX, lnn0_interp, 'r-', 'LineWidth', line_width);  hold on;
	plot(rX, lnn0_fbt, 'b--', 'LineWidth', line_width);
	xline(r_max * rhos0, 'g--', 'LineWidth', line_width);  hold off;
	legend('interp', 'fbt')
	xlabel('r (cm)');
	ylabel('$$\ln n_0/\mathrm{cm^{-3}}$$', 'interpreter', 'latex');
	fig_file = fullfile(fig_path, ['lnn0_', sprintf('%4.4d', idiag), '.png']);
	print(gcf, '-dpng', fig_file);
end

