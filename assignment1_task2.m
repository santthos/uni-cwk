clc
clear
% set initial values and constants 
k = 5.0;
m = 0.5;
T = 10; 
N = 100;
array  = linspace(0, T, N+1);
dt     = T/double(N);
N2 = 1000;
array2  = linspace(0, T, N2+1);
dt2     = T/double(N2);
x0 = 1.0;
v0 = 0.1; 
% call semi-implicit euler method function
u_semi = semi_euler(x0, v0, N, T, k, m); 
u_semi2 = semi_euler(x0, v0, N2, T, k, m); 
% create energy array of zeros
energy = zeros(1,N+1);
% apply for loop to fill energy array for discrete energy 
for i=1:N
    energy(1,i) = 0.5*m*(u_semi(2,i))^2 + 0.5*k*(u_semi(1,i))^2;
end
energy2 = zeros(1,N2+1);
for i=1:N2
    energy2(1,i) = 0.5*m*(u_semi2(2,i))^2 + 0.5*k*(u_semi2(1,i))^2;
end

% apply for loop to fill energy array for discrete modified energy 
energy3 = zeros(1,N+1);
for i=1:N
    energy3(1,i) = 0.5*m*(u_semi(2,i))^2 + 0.5*k*(u_semi(1,i))^2 - (k*dt*u_semi(2,i)*u_semi(1,i));
end
energy4 = zeros(1,N2+1);
for i=1:N2
   energy4(1,i) = (0.5*m*(u_semi2(2,i))^2) + (0.5*k*(u_semi2(1,i))^2) - (k*dt2*u_semi2(2,i)*u_semi2(1,i));
end

% plot figure
figure(1)
plot(array, energy(1,:), 'r-'); hold on;
plot(array2, energy2(1,:), 'b-'); hold on;
plot(array, energy3(1,:), 'g-'); hold on;
plot(array2, energy4(1,:), 'm-'); 
legend('N = 100, Discrete Energy', 'N = 1000, Discrete Energy','N = 100, Discrete Modified Energy','N = 1000, Discrete Modified Energy', 'location', 'NorthWest');
ylim([0 4]);
xlabel('Time','FontSize',11);
ylabel('Energy', 'FontSize', 11);



    