function [t, U_soln] = diffusion3d( kappa, h, U_init, U_bndry, t_int, nt )

% Error Checking
% ==============
%
%   To ensure all values passed in are correct and valid

    %ensure c and h are scalar values
    if ~isscalar(kappa)
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
    %ensure U_init is a 3D matrix
    dim = size(U_init);
    if dim(1)== 1 || dim(2)== 1 || dim(3)== 1
        throw(MException('MATLAB:invalid_argument', ...
        'U_init is not 2D'));
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
    dt = (tf - t0)/(nt - 1);
    t = linspace( t0, tf, nt );
    
    %ensure requirements for convergence are met, give better nt if not.
    r = kappa*dt/h^2;
    if r > 0.25
        nt_star = ceil((kappa*(t_rng(2)-t_rng(1)))/(0.25*h*h) +1);
        throw( MException( 'MATLAB:invalid_argument', ...
        'The ratio kappa*dt/h^2 = %f >= 0.5. Consider using nt = %d', ...
        ratio,nt_star ));
    end

% Solving the system
% ==============
%
%   Use finite difference equation to solve for unknown points.

    [nx, ny, nz] = size( U_init );
 
    %set intial conditions 
    U_soln = zeros( nx, ny, nz, nt );
    U_soln(:, :, :, 1) = U_init;
   
    %iterate through all remaining points in time
    for it = 2:nt
        U_soln(:, :, :, it) = U_bndry( t(it), nx, ny ,nz );
        
        %for each point in the solution matrix
        for ix = 1:nx
            for iy = 1:ny
                for iz = 1:nz
                    %if the point is unknown and is to be solved for
                    if U_soln(ix, iy, iz, it) == -Inf
                        %update current point using points in the past
                        Utmp = U_soln(ix, iy, iz, it - 1);
                        U_soln(ix, iy, iz, it) = Utmp;

                        %define vector to find adjacent points from current
                        %point and save coordinates to dix and diy
                        for dxyz = [-1 1 0 0 0 0; 0 0 -1 1 0 0; 0 0 0 0 -1 1];
                            dix = ix + dxyz(1);
                            diy = iy + dxyz(2);
                            diz = iz + dxyz(3);

                            %if the adjacent point is a finite value then add
                            %term to finite difference equation else do nothing
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
