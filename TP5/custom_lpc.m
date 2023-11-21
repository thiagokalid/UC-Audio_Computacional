function [a, E] = custom_lpc(x, p)
    % Custom LPC function using the autocorrelation method
    N = length(x);
    
    % Calculate autocorrelation function
    R = zeros(1, p+1);
    for k = 0:p
        R(k+1) = sum(x(1:N-k) .* x(1+k:N));
    end
    
    % Toeplitz matrix
    R_matrix = toeplitz(R(1:p));
    
    % Autocorrelation vector
    r = -R(2:p+1)';
    
    % Solve the Yule-Walker equations
    a = linsolve(R_matrix, r);
    
    % Calculate prediction error
    E = R(1) + sum(a .* R(2:p+1));
end