% wave 2d
%
% Parameters
% ==========
% c Some constant in the wave equation
% h the distance between points when approximating the solution or the 
% Delta x between each point in space. 
%
% u_init An nx x ny array array containing initial conditions for the equation/system. 
% u_bndry A nx x ny array of finite values for any Dirichlett value and NaN
% for any insulated boundary, finally –Inf for unknown values to solve for. 
% du_init An nx x ny array array containing the initial rates of change for the wave equation. 
%
% t_int The time interval [t0, tf] on which we will
% approximate the solution.
% nt We will approximate the solution at nt equally spaced time steps for
% t0 <= t <= tf.
%
% Return Values
% =============
% t The vector t contains the nt points in time at which we will be
% approximating the solution.
% U_soln The matrix U is an nx x nt matrix where U(k,t) is the approximation of
% the solution the location x(k) and time t(t).

function [t, U_soln] = wave3d( c, h, U_init, dU_init, U_bndry, t_int, nt )

% Error Checking
% ==============
%
%   To ensure all values passed in are correct and valid

    %ensure c and h are scalar values
    if ~isscalar(c)
        throw(MException('MATLAB:invalid_argument', ...
        'c is not a scalar'));
    end
    if ~isscalar(h) || h < 0
        throw(MException('MATLAB:invalid_argument', ...
        'h is not a positive scalar'));
    end
    %ensure U_bndry is a function handle
    if(~isa( U_bndry, 'function_handle' ))
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument U_bndry is not a function handle' ) );
    end
    %ensure U_init is a 2D matrix
    dim = size(U_init);
    if dim(1)== 1 || dim(2)== 1 || dim(3)== 1
        throw(MException('MATLAB:invalid_argument', ...
        'U_init is not 2D'));
    end
    %ensure dU_init is a 2D matrix
    dim2 = size(dU_init);
    if dim2(1)== 1 || dim2(2)== 1 || dim2(3)== 1
        throw(MException('MATLAB:invalid_argument', ...
        'dU_init is not 2D'));
    end
    %ensure dimentions of U_init and dU_init are the same
    if dim ~= dim2
        throw(MException('MATLAB:invalid_argument', ...
        'U_init, dU_init do not have same dimentions'));
    end
    %ensure t_int is a row vector
    row_dim = size(t_int);
    if row_dim(1) ~= 1 
        throw(MException('MATLAB:invalid_argument', ...
        't_int is not a row vector'));
    end 
    %ensure nt is a positive integer
    if ceil(nt) ~= nt || nt <= 0 || ~isscalar(nt)
        throw(MException('MATLAB:invalid_argument', ...
        'nt must be a positive integer'));
    end
    
    t0 = t_int(1);
    tf = t_int(2);
    t = linspace( t0, tf, nt );
    dt = (tf - t0)/(nt - 1);
    
    %ensure requirements for convergence are met, give better nt if not.
    ratio = (c*dt/h)^2;
    if ratio > 0.25
        nt_star = ceil(((c*(t_rng(2)-t_rng(1)))/(0.5*h)) +1);
        throw( MException( 'MATLAB:invalid_argument', ...
        'The ratio c*dt/h^2 = %f >= 0.5. Consider using nt = %d', ...
        ratio,nt_star ));
    end
 
% Initialization
% ==============
%
%   Your description here.

    [nx, ny, nz] = size( U_init );
 
    U_soln = zeros( nx, ny, nz, nt );
    U_soln(:, :, :, 1) = U_init;
    
    r = (c*dt/h)^2;

    U_soln(:, :, :, 2) = U_soln(:, :, :, 1) + dt*dU_init;

    for it = 3:nt
        U_soln(:, :, :, it) = U_bndry( t(it), nx, ny, nz);
        
        for ix = 1:nx
            for iy = 1:ny
                for iz = 1:nz
                    if U_soln(ix, iy, iz, it) == -Inf
                        Utmp = U_soln(ix, iy, iz, it - 1);
                        U_soln(ix, iy, iz, it) = 2*Utmp - U_soln(ix, iy, iz, it - 2);

                        for dxyz = [-1 1 0 0 0 0; 0 0 -1 1 0 0; 0 0 0 0 -1 1];
                            dix = ix + dxyz(1);
                            diy = iy + dxyz(2);
                            diz = iz + dxyz(3);

                            if ~isnan( U_soln(dix, diy, diz, it - 1) )
                                U_soln(ix, iy, iz, it) = U_soln(ix, iy, iz, it) + ...
                                    r*( U_soln(dix, diy, diz, it - 1) - Utmp );
                            end
                        end
                    end
                end
            end
        end
    end
end

