clear all;
close all;

tend = 1.0;
nsteps = 10;

v0 = 1.0; % set transport velocity

% Number of mesh points
N = 64;
L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);
dx    = xmesh(2) - xmesh(1);

% create finite difference matrix
e = ones(N,1);
D = spdiags([-e 0*e e], -1:1, N, N);
D(1,N) = -1.0;
D(N,1) = 1.0;

D = 1.0/(2.0*dx)*D;

% define initial value
u_0  = @(x) sin(x);
u0 = u_0(xmesh);


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
    plot(xmesh, u(:,n), 'b'); hold on;
    plot(xmesh, u_0(xmesh), 'r');
    xlim([xmesh(1), xmesh(end)]);
    ylim([-1.1, 1.1]);
    xlabel('x');
    ylabel('u');
    title('Solution in physical space');
    txt = strcat('t=',num2str(taxis(n), '%3.2f'));
    text(0.75*L, 0.75, txt);
    drawnow;
    
%     if n==length(taxis)
%        filename = 'transport_solution.pdf';
%        export_fig(filename);
%     end
    
end

figure(2);
plot(real(eig(-v0*full(D))), imag(-v0*eig(full(D))), 'bo', 'markerfacecolor', 'b');
xlim([-1.0, 1.0]);
ylim([-20, 20]);
xlabel('Real part of eigenvalue');
ylabel('Imaginary part of eigenvalue');