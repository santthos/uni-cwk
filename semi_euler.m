% function for semi-implicit euler method that takes the initial values x0,
% v0, the number of time steps N, the final time T and the parameters k and
% m as arguments

function u = semi_euler(x0, v0, nsteps, tend, k, m)
dt = tend/double(nsteps);
u = zeros(2,nsteps+1);
u0 = [x0;v0];
u(:,1) = u0;;
for i=1:nsteps
    u(2,i+1) = u(2,i) - ((u(1,i))*dt*(k/m));
    u(1,i+1) = u(1,i) + (u(2,i)*dt) - ((dt)^2*(k/m)*u(1,i));   
end
end  