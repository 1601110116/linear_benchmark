% For a given case, the basis functions of Fourier-Bessel series
%   are fixed. This script generate them for further use

% Up to the order of m_max
m_max = 20;
% First l_max positive zeros
l_max = 30;
% The output file
output = 'fb_cores.mat';

bessel_cores = cell(m_max+1, l_max);
fourier_cores = cell(m_max+1, 1);
load('parameters.mat');
addpath(code_path);
build_grid_2d;
for m = 0: m_max
	disp(['m = ', num2str(m), ' of ', num2str(m_max)]);
	for l = 1: l_max
		lambda_ml = alpha(m+1, l) / radius;
		bessel_cores{m+1, l} = besselj(m, lambda_ml * r2d);
	end
	fourier_cores{m+1} = exp(-1j * m * tht2d);
end
save(output, 'bessel_cores', 'fourier_cores');
