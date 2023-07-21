clear all;
close all;

%tend = 3.55; % try values of 2.0 (before shock formation) and 4.0 (after shock formation)
tend = 3.0;
%nu   = 0.002; % try values 0.002 (weak diffusion) and 1.0 (strong diffusion)
nu   = 0.002;

% Set number of Fourier modes and generate mesh; as always, throw away last
% point which is identical to the first one because of the assumed periodic
% boundary conditions
N = 64;
assert(mod(N,2)==0, 'N must be even');
L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);

% create spectral derivative matrix 
D = diag(1j*[0:N/2 -N/2+1:-1]);

% define initial value
u_0  = @(x) sin(x);

% compute u by evaluating u0 on the mesh and transforming into spectral
% space
u0 = fft(u_0(xmesh).');

% define the RHS. We use Matlab's fft and ifft function to transfer back
% and forth between physical and spectral space.
% Note how we resolve the nonlinearity by first computing the derivative in
% spectral space through D_x*u, transform u and u_x into physical space using ifft,
% computing u.*u_x and then transforming back into spectral space using
% fft. The linear diffusion term is computed completely in spectral space
% and requires no fft/ifft.
f = @(t,u) -fft( ifft(u).*ifft(D*u) ) + nu*D*D*u;
sol = ode45(f, [0 tend], u0);

fprintf('Minimum time step used: %5.3e \n', min(diff(sol.x)));

% for plotting, interpolate solution to a time mesh we define
taxis = linspace(0, tend, 200);
u = deval(taxis, sol);

% print out solutions at t=0
figure(1)
clf
subplot(211), plot(xmesh, real(ifft(u(:,1))), 'r'); hold on;
xlim([xmesh(1), xmesh(end)]);
ylim([-1.1, 1.1]);
xlabel('x');
ylabel('u');
title('Solution in physical space');
txt = strcat('t=',num2str(0.0, '%3.2f'));
text(0.75*L, 0.75, txt);

% plot solution in spectral space; we use fftshift to move the constant
% wave number mode (k=0) into the centre.
subplot(212), plot(-N/2:N/2-1, abs(fftshift(u(:,1)))/N, 'rd', 'markerfacecolor', 'r');
xlim([-N/2, N/2-1]);
ylim([0.0, 0.55]);
xlabel('k');
ylabel('u_k');
title('Solution in spectral space');

% print out solution at all times
figure(2);
for n=1:length(taxis)
    clf;
    
    % plot solution in physical space by applying ifft and taking the real
    % part - note how the solution becomes very steep for small viscosity
    % until a discontinuity in the first derivative forms (a "shock")
    subplot(211), plot(xmesh, real(ifft(u(:,n))), 'b'); hold on;
    plot(xmesh, real(ifft(u(:,1))), 'r');
    xlim([xmesh(1), xmesh(end)]);
    ylim([-1.1, 1.1]);
    xlabel('x');
    ylabel('u');
    title('Solution in physical space');
    txt = strcat('t=',num2str(taxis(n), '%3.2f'));
    text(0.75*L, 0.75, txt);
    
    % plot solution in spectral space; we use fftshift to move the constant
    % wave number mode (k=0) into the centre.
    subplot(212), plot(-N/2:N/2-1, abs(fftshift(u(:,1)))/N, 'rd', 'markerfacecolor', 'r');
    plot(-N/2:N/2-1, abs(fftshift(u(:,n)))/N, 'bo', 'markerfacecolor', 'b');
    xlim([-N/2, N/2-1]);
    ylim([0.0, 0.55]);
    xlabel('k');
    ylabel('u_k');
    title('Solution in spectral space');
    drawnow;
end