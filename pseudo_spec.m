%input variable u_0 as 'f= @(x) ...' 
%the function outputs the solution in physical space which can be plotted
%using plot(xmesh, u, 'b') after running the function
function [u, xmesh] = pseudo_spec(N, tend, v_0, nu, u_0)

%creates the x axis matrix
nsteps = 264;
L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);
%creates the differentiation matrix
D = diag(1j*[0:N/2 -N/2+1:-1]);

u0 = fft(u_0(xmesh).');
%the function as written compactly with the differentiation matrix
f = @(t,u) (-(v_0*D*u)+(nu*D*D*u));
taxis = linspace(0, tend, nsteps+1);

%this is the explicit euler solver function with 4x the number of fourier
%modes as the time steps
nsteps2 = N*4 ;
dt   = tend/double(nsteps2);
u2    = zeros(length(u0),nsteps2+1);
u2(:,1) = u0;
fcount = 0;
for i=1:nsteps
    u2(:,i+1) = u2(:,i) + dt*f(i*dt, u2(:,i));
    fcount = fcount + 1;
end
%prints solutions in physical space 
u = (real(ifft(u2(:,nsteps+1))));

