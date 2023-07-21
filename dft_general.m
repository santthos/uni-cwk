clear all;
close all;

% number of discrete Fourier modes
Nx    = 25;
xmesh = linspace(0, 2*pi, Nx+1);
% remove right boundary point
xmesh = xmesh(1:Nx);

u = @(x) heaviside(x-2*pi/3) - heaviside(x - 4*pi/3);

% Compute FFT on normal mesh
uhat = fft(u(xmesh));

% Transform back to physical space
uifft = ifft(uhat);

% Create a fine mesh to visualise trunaced part of spectrum
fact    = 40; % resolution factor between coarse and fine mesh
Nx_fine = fact*Nx + 1;
xmesh_fine = linspace(0, 2*pi, Nx_fine+1);
% remove last point b/c of periodicity
xmesh_fine = xmesh_fine(1:end-1);
% find u on fine mesh and DFT on fine mesh
u_fine     = u(xmesh_fine);
u_fine_hat = fft(u_fine);

% Compute the inverste FFT to provide values on the fine mesh; need to
% normalise the amplitudes to correct for the different number of values
% Nx vs. Nx_fine
ureconstructed = (Nx_fine/Nx)*ifft(uhat, Nx_fine);

figure(1);
plot(xmesh_fine, real(ureconstructed), 'r'); hold on;
plot(xmesh_fine, real(u_fine), 'b--');
plot(xmesh,      real(uifft), 'ro', 'markerfacecolor', 'r');
legend('DFT representation','Actual function');
title('Functions in physical space');
xlabel('x');
ylabel('u');

