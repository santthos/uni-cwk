% 
% This script generates the matrix D that represents our finite difference
% approximation of the double spatial derivative.
%
% We then test it on a simple function and check it does actually
% approximate the second derivative.
%

clear all;
close all;

L = 2.3;
nnodes = [10, 25, 50, 75, 100, 250, 500, 1000];
nus    = [0.1, 1.0, 5.0];
eigvals = zeros(3, length(nnodes));
dxs    = zeros(3, length(nnodes));

for m=1:length(nus)
    nu = nus(m);
    for n=1:length(nnodes)
        N = nnodes(n);
        xaxis = linspace(0, L, N+1); % create nodes
        dx = xaxis(2) - xaxis(1);
        dxs(m, n) = dx;
        % Create the matrix D using spdiags
        e = ones(N-1,1);
        D = (nu/dx^2)*spdiags([e -2*e e], -1:1, N-1, N-1);
        
        % Store largest eigenvalue of D
        eigvals(m,n) = max(abs(eig(D)));
    end
end

figure(1);
loglog(nnodes, eigvals(1,:), 'bo', 'markerfacecolor', 'b'); hold on;
loglog(nnodes, eigvals(2,:), 'rs', 'markerfacecolor', 'r');
loglog(nnodes, eigvals(3,:), 'gp', 'markerfacecolor', 'g');
% The 4 is a constant and only approximate ... gives a very good fit though
loglog(nnodes, 4*nus(1)./dxs(1,:).^2, 'b-');
loglog(nnodes, 4*nus(2)./dxs(2,:).^2, 'r-');
loglog(nnodes, 4*nus(3)./dxs(3,:).^2, 'g-');
legend(strcat('\nu = ', num2str(nus(1), '%3.2f')), strcat('\nu = ', num2str(nus(2), '%3.2f')), strcat('\nu = ', num2str(nus(3), '%3.2f')), 'location', 'NorthWest');
xlabel('Number of nodes N');
ylabel('Largest eigenvalue of \nu D');
