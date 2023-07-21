function [u, fcount] = exp_euler(u0, tend, nsteps, f)
dt   = tend/double(nsteps);
u    = zeros(length(u0),nsteps+1);
u(:,1) = u0;
fcount = 0;
for i=1:nsteps
    u(:,i+1) = u(:,i) + dt*f(i*dt, u(:,i));
    fcount = fcount + 1;
end
end