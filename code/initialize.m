function initialize (simulate_mode, N0, Ln, init_perturbation)
% This function gets the inittial values of all simulated quantities
global Te den Te_aux den_aux vi w vi_aux w_aux jz ve x ...
	phi vEx vEy last_diagnose nx nz calc vdex vdey data_path
jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz;  phi = jz;

if simulate_mode == 1
	x3d1 = repmat(reshape(x, [], 1, 1), 1, nx+2, nz+2);
	y3d1 = repmat(reshape(x, 1, [], 1), nx+2, 1, nz+2);
	r3d1 = sqrt(x3d1.^2 + y3d1.^2);

	Te = ones(nx+2, nx+2, nz+2);
	den = N0 * exp(-0.5 * (r3d1 / Ln) .^ 2);
%	Te(2:end-1, 2:end-1, 2:end-1) = Te(2:end-1, 2:end-1, 2:end-1) + ...
%		init_perturbation * calc .* (rand(nx, nx, nz) - 0.5);
%	den(2:end-1, 2:end-1, 2:end-1) = den(2:end-1, 2:end-1, 2:end-1) + ...
%		init_perturbation * calc .* (rand(nx, nx, nz) - 0.5);
	Te(2:end-1, 2:end-1, 2:end-1) = Te(2:end-1, 2:end-1, 2:end-1) + ...
		init_perturbation * calc .* (rand(nx, nx, nz));
	den(2:end-1, 2:end-1, 2:end-1) = den(2:end-1, 2:end-1, 2:end-1) + ...
		init_perturbation * calc .* (rand(nx, nx, nz));
	w = zeros(nx+2, nx+2, nz+2);  vi = w;
%	w(2:end-1, 2:end-1, 2:end-1) = init_perturbation * calc .* (rand(nx,nx,nz)-0.5);
	sphi(nx, nz);
	sdata();
elseif simulate_mode == 2
	last_file = get_last_file(data_path);
	disp(['started from file ', last_file]);
	load(last_file)
	sdata();
elseif simulate_mode == 3
	last_file = get_last_file(fullfile('..', 'data'));
	disp(['started from ', last_file]);
	last_diagnose = 0;
	% copy the last data file in the previous directory to current
	% directory as rest.mat and then load it
	copyfile(last_file, 'rest.mat');
	load rest.mat
	sdata();
end
den_aux = den;  Te_aux = Te;
vi_aux = vi;      w_aux = w;

