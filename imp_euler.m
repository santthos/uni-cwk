function u = imp_euler(u0, tend, nsteps, lambda)
dt   = tend/double(nsteps);
u    = zeros(1,nsteps+1);
u(1) = u0;
for i=1:nsteps
    u(i+1) = u(i)/(1 + dt*lambda);
end
end