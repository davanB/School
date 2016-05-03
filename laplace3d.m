% laplace3D
% In 2 dimensions our finite difference equations for approximating
% Laplaces equation refers to points on grid that can be used to build a
% matrix and calculate the interior points. In 3 dimensions this reflects 
% this turns into a volume such as a cube. The boundaries of the surfaces
% or shapes correspond to the boundary values of the problem. 
%
% Parameters
% ==========
%    U - The 3D matrix specifying the boundary conditions 
%        and which points to approximate the solution.
%
% Return Values
% =============
%    U_out - The matrix Uout is a 3D matrix which contains the 
%            approximation o the solution to Laplace’s equation.

function [U_out] = laplace3d( U )
    [n_x, n_y, n_z] = size( U );
    U_out = U;
    
    %Error checking 
    %ensure all boundry points are not -Inf
    if n_x < 2 || n_y < 2 || n_z < 2
        throw( MException( 'MATLAB:invalid_argument', ...
                'the matrix passed in is not 3D' ) );
    end
    for z = 1:n_z
        for n = 1:n_x
            if U(n,1,z) == -Inf || U(n,end,z) == -Inf
                throw( MException( 'MATLAB:invalid_argument', ...
                'the matrix has a boundry value of -Inf' ) );
            end
        end
        for j = 1:n_y
            if U(1,j,z) == -Inf || U(end,j,z) == -Inf
                 throw( MException( 'MATLAB:invalid_argument', ...
                'the matrix has a boundry value of -Inf' ) );
            end
        end
    end
    
    % Step 2
    u_to_w = zeros( n_x, n_y, n_z );
    w_to_u = zeros( 3, n_x * n_y * n_z );
    m = 0;

    for ix = 1:n_x
        for iy = 1:n_y
            for iz = 1:n_z
                if U(ix, iy, iz) == -Inf
                    m = m + 1;
                    u_to_w(ix, iy, iz) = m;
                    w_to_u(:, m) = [ix, iy, iz]';
                end
            end
        end
    end

    % Create the sparse system of linear equations
    M = spalloc( m, m, 7*m );
    b = zeros( m, 1 );

    for k = 1:m
        % Get the coordinates of the kth point
        % For each of the 6 adjacent points, determine if
        % the point is an insluated boundary point, a Dirichlet
        % boundary point or an unknown value and modify M as appropriate.
        c = w_to_u(:,k);
        c1 = c + [-1  0  0]';
        p1 = U(c1(1),c1(2),c1(3));
        c2 = c + [ 1  0  0]';
        p2 = U(c2(1),c2(2),c2(3));
        c3 = c + [ 0 -1  0]';
        p3 = U(c3(1),c3(2),c3(3));
        c4 = c + [ 0  1  0]';
        p4 = U(c4(1),c4(2),c4(3));
        c5 = c + [ 0  0 -1]';
        p5 = U(c5(1),c5(2),c5(3));
        c6 = c + [ 0  0  1]';
        p6 = U(c6(1),c6(2),c6(3));
        
        if p1 ~= -Inf && ~isnan(p1)%Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p1;
        elseif p1 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(c1(1),c1(2),c1(3));
            M(k,j) = M(k,j) + 1;
        end
        if p2 ~= -Inf && ~isnan(p2) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p2;
        elseif p2 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(c2(1),c2(2),c2(3));
            M(k,j) = M(k,j) + 1;
        end    
        if p3 ~= -Inf && ~isnan(p3) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p3;
        elseif p3 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(c3(1),c3(2),c3(3));
            M(k,j) = M(k,j) + 1;
        end
        if p4 ~= -Inf && ~isnan(p4) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p4;
        elseif p4 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(c4(1),c4(2),c4(3));
            M(k,j) = M(k,j) + 1;
        end
         if p5 ~= -Inf && ~isnan(p5) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p5;
        elseif p5 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(c5(1),c5(2),c5(3));
            M(k,j) = M(k,j) + 1;
         end
         if p6 ~= -Inf && ~isnan(p6) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p6;
        elseif p6 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(c6(1),c6(2),c6(3));
            M(k,j) = M(k,j) + 1;
        end
    end

    w = M \ b;

    for k = 1:m
        % Copy the value from w into U_out
        point = w_to_u(:,k);
        value = w(k);
        U_out(point(1),point(2),point(3)) = value;
    end
end
