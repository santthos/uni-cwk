clc
clear
% set initial values and constants 
k = 5.0;
m = 0.5;
T = 10; 
N = 100;

x0 = 1.0;
v0 = 0.1; 
u0 = [x0;v0];

array  = linspace(0, T, N+1);

f = @(t,u) [u(2); (-(k/m)*u(1))]; 

[time_ode45, u_ode45] = ode45(f, [0, T], u0);

%call semi-implicit euler method function
u_semi = semi_euler(x0, v0, N, T, k, m);

% creat figure 
set(gcf,'Units','centimeter',  'Position',[0 0 18 6]);

set(gcf, 'Color', 'White');

figure(1);
plot(array, u_semi(1,:), 'r', 'LineWidth', 2); hold on;
plot(time_ode45, u_ode45(:,1), 'k', 'LineWidth', 2);
legend('semi-implicit', 'ode45', 'Location','SouthWest');
ylim([-2.5, 2.5]);
xlim([array(1) array(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);



    