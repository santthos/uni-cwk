clear all;
close all;

% set final time and diffusivity
tend = 1.0;
nu   = 0.5;

% number of Fourier modes to use
Nx    = 48;
assert(mod(Nx,2)==0, 'Nx must be even');

% As always, we use [0, 2*pi] as domain and throw away the very last mesh
% point because it is identical to the first in periodic cases
L = 2*pi;
mesh = linspace(0, L, Nx+1);
mesh = mesh(1:end-1);

% create derivative matrix - note that we use the "right" matrix where
% instead of going from 0:N-1, we go from 0 to N/2-1 and then start from
% negative -N/2+1 and count upwards to -1 to avoid issues with computing
% derivatives of aliased modes
% create spectral derivative matrix 
D = diag(1j*[0:Nx/2 -Nx/2+1:-1]);
Dsquared = diag(-[0:Nx/2-1 0 -Nx/2+1:-1].^2);

% define initial value
n   = 2;
lam = (n*pi/L)^2;
Tex = @(t,x) exp(-nu*lam*t)*sin(sqrt(lam)*x);

% get mesh values in physical space
T0  = @(x) Tex(0,x);

% compute phi_0 by using fft to transform into spectral space
That_0 = fft(T0(mesh).');

% RHS function in spectral space
f = @(t,x) nu*D*D*x;
[t, y] = ode45(f, [0 tend], That_0);

fprintf('Smallest time step used in simulation: %5.3e \n', min(diff(t)));

figure(1);
for n=1:10:length(t)
    clf;
    T_hat = y(n,:).';
    plot(mesh, real(ifft(T_hat)), 'b'); hold on;
    plot(mesh, Tex(t(n), mesh), 'k+');
    xlim([mesh(1), mesh(end)]);
    ylim([-1.1, 1.1]);
    drawnow;
end