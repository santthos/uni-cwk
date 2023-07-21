clear all;
close all;

tend = 2*pi;
tend = 3.0;
nsteps = 100;

v0 = 1.0; % set transport velocity

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

f = @(t,u) -v0*D*u;
taxis = linspace(0, tend, nsteps+1);

solver = 'exp_euler';

if strcmp(solver, 'ode45')
    
    sol = ode45(f, [0 tend], u0);
    
    fprintf('Minimum time step used: %5.3e \n', min(diff(sol.x)));
    
    % for plotting, interpolate solution to a time mesh we define
    u = deval(taxis, sol);
    
elseif strcmp(solver,'imp_euler')
    % use implicit Euler instead of ode45
    u = imp_euler_linear(u0, tend, nsteps, -v0*D);
    
elseif strcmp(solver,'exp_euler')
    u = exp_euler(u0, tend, nsteps, f);
end

figure(1);
for n=1:length(taxis)
    
    % plot solution in physical space by applying ifft and taking the real
    % part - note how the solution becomes very steep for small viscosity
    % until a discontinuity in the first derivative forms (a "shock")
    figure(1);
    clf;
    plot(xmesh, real(ifft(u(:,n))), 'b'); hold on;
    plot(xmesh, u_0(xmesh), 'r');
    xlim([xmesh(1), xmesh(end)]);
    ylim([-1.1, 1.1]);
    xlabel('x');
    ylabel('u');
    title('Solution in physical space');
    txt = strcat('t=',num2str(taxis(n), '%3.2f'));
    text(0.75*L, 0.75, txt);
    
%     if n==length(taxis)
%        filename = 'transport_solution.pdf';
%        export_fig(filename);
%     end
    
    % plot solution in spectral space; we use fftshift to move the constant
    % wave number mode (k=0) into the centre.
    figure(2);
    clf;
    plot(-N/2:N/2-1, abs(fftshift(u(:,n)))/N, 'bo', 'markerfacecolor', 'b'); hold on;
    plot(-N/2:N/2-1, abs(fftshift(u0))/N, 'rd', 'markerfacecolor', 'r');
    xlim([-N/2, N/2-1]);
    ylim([0.0, 0.55]);
    xlabel('k');
    ylabel('u_k');
    title('Solution in spectral space');
    drawnow;
    
%     if n==length(taxis)
%         filename = 'transport_solution_spectrum.pdf';
%         export_fig(filename);
%     end
end

figure(3);
plot(real(eig(-v0*D)), imag(eig(-v0*D)), 'bo', 'markerfacecolor', 'b');
xlabel('Real part of eigenvalue');
ylabel('Imaginary part of eigenvalue');