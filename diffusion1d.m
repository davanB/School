function [x_out, t_out, U_out] = diffusion1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )

% Error Checking
% ==============
%
%   To ensure all values passed in are correct and valid

%check to make sure x range is a row vector and b > a
x_size = size(x_rng);
if(x_size(1) ~= 1 || x_rng(2) <= x_rng(1))
    throw(MException('MATLAB:invalid_argument', ...
        'invalid x_rng row vector'));
end
%check to make sure t range is a row vector and t_final > t_o
t_size = size(t_rng);
if(t_size(1) ~= 1 || t_rng(2) <= t_rng(1))
    throw(MException('MATLAB:invalid_argument', ...
        'invalid t_rng row vector'));
end
%check to make sure nx is a integer and is greater than 0
if(~isscalar(nx) || nx ~= round(nx) || nx <= 0)
    throw(MException('MATLAB:invalid_argument', ...
        'nx is either not a scalar value or is not an integer greater than zero'));
end
%check to make sure nt is a integer and is greater than 0
if(~isscalar(nt) || nt ~= round(nt) || nt <= 0)
    throw(MException('MATLAB:invalid_argument', ...
        'nt is either not a scalar value or is not an integer greater than zero'));
end
%check to make sure u_init is a function handle
if(~isa( u_init, 'function_handle' ))
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_init is not a function handle' ) );
end
%check to make sure _bndry is a function handle
if(~isa( u_bndry, 'function_handle' ))
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a function handle' ) );
end
%check to ensure the ratio is not greate than 0.5 to ensure proper
%heat diffision equation behavior
delta_t= (t_rng(2) - t_rng(1))/(nt-1);
h = (x_rng(2)-x_rng(1))/(nx-1);
ratio = (kappa*delta_t)/(h*h);
if(ratio > 0.5)
        nt_star = ceil((kappa*(t_rng(2)-t_rng(1)))/(0.5*h*h) +1);
        throw( MException( 'MATLAB:invalid_argument', ...
        'The ratio kappa*dt/h^2 = %f >= 0.5. Consider using nt = %d', ...
        ratio,nt_star ));
end

% Initialization
% ==============
%
%   Begin setup for solving the diffusion equation.

%create a column vector of x values and a row vector of t values
x_column = linspace(x_rng(1),x_rng(2),nx)';
t_row = linspace(t_rng(1), t_rng(2), nt);

%create an empty matrix
matrix = zeros(nx,nt);

%populate the first colunm with initial u_init(xn)
first_col = u_init(x_column);
matrix(:,1) = first_col;

%populate first and last row with b_bndry(tn)
boundries = u_bndry(t_row(2:end));
first_row = boundries(1);
last_row = boundries(2);
matrix(1,2:nt) = first_row;
matrix(nx,2:nt) = last_row;

% Solving
% =======
%
%   Iterate through the matrix and polulate fields to solve equation.

%formula to populate matrix with u(i,k) values
for k = 1:(nt-1)
    matrix(2:(nx-1), k+1) = matrix(2:(nx-1),k) + ratio*diff(matrix(:,k),2);
end

x_out = x_column;
t_out = t_row;
U_out = matrix;

end

