function [U, H, its] = poldec(A)
%POLDEC     Polar Decomposition
%   [U, H, ITS] = poldec(A) computes the polar decomposition A = U*H
%   of the square, nonsingular matrix A. ITS is the number of 
%   iterations for convergence.

plot_display = true; %Set to true to plot relative error vs iterations

if plot_display
    errors = [1];
end

X = A;
X_old = zeros(size(X));
tol = sqrt(size(X))*eps;    %Tolerance to test convergence

I_size = eye(size(X));  %Identity matrix of appropriate size
its = 0;
theta = 0.7; %parameter to switch to NS

ns_switch = false;  %Boolean to determine switch to NS

delta = 1;
delta_old = inf;
X_product = X'*X;
residual_inf = norm(I_size - X_product, inf);

cond = false; %if true will test delta < delta_old/2 for convergence

while 1
    its = its + 1;
    fprintf('---\nIteration: %2.0f\n', its)
    X_old = X;    
    
    if ~ns_switch
        if residual_inf < theta
            ns_switch = true;
        end
    end    
    
    if ns_switch %If NS will converge do, otherwise Newton        
        X = (3/2)*X - X*X_product/2;  %Newton-Schulz
        fprintf("Switched to NS\n");
        X_product = X'*X; %MATLAB should detect symm in this
    else
        X = (X + (inv(X))')/2;    %Newton
        X_product = X'*X;
    end            
    
    delta_old = delta;
    delta = norm(X - X_old, inf)/norm(X, inf);
    fprintf('Relative Error: %g\n', delta)
    if plot_display
        errors(its) = delta;
    end
    residual_inf = norm(I_size - X_product, inf);
    fprintf('Residual: %g\n', residual_inf)
    
    if delta < 1e-3 %Don't want early termination when delta near 1
        cond = true;
    end
    
    if (delta < tol) | ((delta > delta_old/2) & cond)
        break
    end
    
end

if plot_display
    semilogy(1:its, errors, '-o')
    xlabel('Iterations')
    leg = legend('$\delta_k$')
    set(leg, 'Interpreter', 'Latex')
    set(leg, 'FontSize', 12)
    set(leg, 'Location', 'best')
end

U = X;
H = U'*A;
H = 0.5*(H + H'); %H computed naively not Hermitian in general


