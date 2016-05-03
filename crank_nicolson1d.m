% crank_nicolson1d
%
% We are solving the heat-conduction/diffusion equation on an interval
% a <= x <= b for a time interval t0 <= t <= tf. The initial conditions are
% given by the function u_init:[a, b] -> R and the boundary conditions are
% given by two functions a_bndry:[t0, t_f] -> R and
% b_bndry:[t0, tf] -> R.
%
% The solution U(x, t) must satisfy the initial and boundary conditions and
% must satisfy the heat-conduction/diffusion equation at every point in the
% region. We will approximate the solution on an nx x nt grid of points.
%
% Parameters
% ==========
% kappa The diffusion coefficient dictated by the system.
% x_rng The spacial interval [a, b] on which we will approximate the
% solution.
% t_rng The time interval [t0, tf] on which we will
% approximate the solution.
%
% u_init A real-valued function of a real variable describing the
% initial condition of the property at time t0nitial.
% u_bndry A vector-valued function returning a 2-dimensional column vector
% of the boundary values a_bndry(t) and b_bndry(t).
%
% nx We will approximate the solution at nx equally spaced points on [a, b].
% nt We will approximate the solution at nt equally spaced time steps for
% t0 <= t <= tf.
%
% Return Values
% =============
% x_out The vector x contains the nx points at which we are approximating the
% solution in the spatial dimension.
% t_out The vector t contains the nt points in time at which we will be
% approximating the solution.
% U_out The matrix U is an nx x nt matrix where U(k,ell) is the approximation of
% the solution the location x(k) and time t(ell). 

function [x_out, t_out, U_out] = ...
       crank_nicolson1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )


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
if(~isscalar(nx) || nx ~= round(nx) || nx <= 4)
    throw(MException('MATLAB:invalid_argument', ...
        'nx is either not a scalar value or is not an integer greater than four'));
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
%warn user if ratio is above 0.5 
delta_t= (t_rng(2) - t_rng(1))/(nt-1);
h = (x_rng(2)-x_rng(1))/(nx-1);
ratio = (kappa*delta_t)/(h*h);
if(ratio > 0.5)
    warning( 'MATLAB:questionable_argument', ...
    'the arguments of %d and %d are sub-optimal', ...
     x_rng(1), x_rng(2) )
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

%populate first and last row with u_bndry(tn)
for i = 2:nt
    temp_mat = u_bndry(t_row(i));
    matrix(1,i) = temp_mat(1);
    matrix(nx,i) = temp_mat(2);
end

% Solving
% =======
%
%   Iterate through the matrix and polulate fields to solve equation.

%formula to populate matrix with u(i,k) values
for k = 1:(nt-1)
    %set up tri diagonal matrix 
    Tridiagonal_matrix1 = diag( -1*ratio * ones(nx-3, 1), 1);   
    Tridiagonal_matrix2 = diag( -1*ratio * ones(nx-3,1), -1);
    Tridiagonal_matrix3 = diag( 2*(1+ratio)*ones(nx-2, 1) );
    if isnan(matrix(1,k+1)) == 1
        Tridiagonal_matrix2(1,2) = (-2/3)*ratio;
        Tridiagonal_matrix3(1,1) = 2+2/3*ratio;
    end
    if isnan(matrix(nx,k+1)) == 1
        Tridiagonal_matrix2(end,end-1) = -2/3*ratio;
        Tridiagonal_matrix3(end,end) = 2+2/3*ratio;
    end
    %complete tri-diagonal matrix 
    Tridiagonal_matrix = Tridiagonal_matrix1 + Tridiagonal_matrix2 + Tridiagonal_matrix3;

    %boundries cond is simply to easily add to the "b" column vector in
    %Mx=b
    boundries_cond = zeros(nx-2,1);
    boundries_cond(1) = ratio*matrix(1,k+1);
    boundries_cond(end) = ratio*matrix(end,k+1);
    if isnan(matrix(1,k+1)) == 1
        boundries_cond(1) = 0;
    end
    if isnan(matrix(nx,k+1)) == 1
        boundries_cond(end) = 0;
    end
    % b is 2 times the value in current column + ratio*the difference eqn
    % the fisrt and last row also get the boundry conditions added to them
    b_known = 2.*matrix(2:(nx-1),k) + ratio.*diff(matrix(:,k),2) + boundries_cond;
    matrix(2:(nx-1),k+1) = (Tridiagonal_matrix\b_known);
    if isnan(matrix(1,k+1)) == 1
        matrix(1,k+1) = 4/3*matrix(2,k+1) - 1/3*matrix(3,k+1);
    end
    if isnan(matrix(nx,k+1)) == 1
        matrix(nx,k+1) = 4/3*matrix(nx-1,k+1) - 1/3*matrix(nx-2,k+1);
    end
end

x_out = x_column;
t_out = t_row;
U_out = matrix;

end

