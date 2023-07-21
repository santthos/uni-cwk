% function for semi-implicit euler method for the nonlinear system that takes the initial values x0,
% v0, the number of time steps N, the final time T, the parameters k and
% m and beta as arguments

function u = semi_euler_nonlinear(x0, v0, nsteps, tend, k, m, beta)
dt = tend/double(nsteps);
u = zeros(2,nsteps+1);
u0 = [x0;v0];
u(:,1) = u0;
for i=1:nsteps
    u(2,i+1) = u(2,i) - ((u(1,i)+((u(1,i))^3*beta))*dt*(k/m));
    u(1,i+1) = u(1,i) + (u(2,i)*dt) - ((dt)^2*(k/m)*u(1,i)) -((dt)^2*((k*beta)/m)*(u(1,i))^3) ;   
end
end  