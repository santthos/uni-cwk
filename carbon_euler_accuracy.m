clc
clear

% Set a few parameter
T      = 1.0;         % final time until which we compute
lambda = 1.0;        % set the decay constant

r0    = 1.0; % we set the ratio at time t=0 to 1.0

% define the exact solution as a function handle
u_exact = @(t) r0*exp(-lambda*t);

% it is easier to define the number of time steps than the dt directly
N = [1000, 750, 500, 250, 100, 75, 50, 10];

% allocate vectors to store errors for every run
err_exp = zeros(1,length(N));
err_imp = zeros(1,length(N));
dts     = zeros(1,length(N));

for n=1:length(N)
    taxis = linspace(0, T, N(n)+1); % add +1 to account for t=0
    u_exp = exp_euler(r0, T, N(n), lambda);
    u_imp = imp_euler(r0, T, N(n), lambda);
    
    % store the time step dt for plotting
    dts(n) = taxis(2) - taxis(1);
    
    % Now compute the errors
    err_exp(1,n) = max(abs(u_exp - u_exact(taxis)));
    err_imp(1,n) = max(abs(u_imp - u_exact(taxis)));
end

% We fit a line log(err) = p*log(N) + C through the data points for
% reasons that will become clear later.
p_exp = polyfit(log(dts), log(err_exp), 1);
p_imp = polyfit(log(dts), log(err_imp), 1);

figure(1);
loglog(dts, err_exp, 'ro', 'markerfacecolor', 'r'); hold on;
% Because we fitted log(err), we need to apply exp to obtain err.
loglog(dts, exp(polyval(p_exp, log(dts))), 'r-');
txt = strcat('Slope p=', num2str(p_exp(1),'%3.2f'));
text(1e-2, 1e-3, txt); % this needs retuning if values in N are changed
xlim([dts(1) dts(end)])
xlabel('\Delta t');
ylabel('Error');
grid on;
legend('Explicit Euler', 'Linear fit', 'location', 'NorthWest');

figure(2)
loglog(dts, err_imp, 'bs', 'markerfacecolor', 'b'); hold on;
loglog(dts, exp(polyval(p_imp, log(dts))), 'b-');
txt = strcat('Slope p=', num2str(p_imp(1),'%3.2f'));
text(1e-2, 1e-3, txt); % this needs retuning if values in N are changed
xlim([dts(1) dts(end)])
xlabel('\Delta t');
ylabel('Error');
grid on;
legend('Implicit Euler', 'Linear fit', 'location', 'NorthWest');

