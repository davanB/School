% laplace2d
% In 2 dimensions our finite difference equations for approximating
% Laplaces equation refers to points on grid that can be used to build a
% matrix and calculate the interior points. In 3 dimensions this reflects 
% this turns into a volume such as a cube. The boundaries of the surfaces
% or shapes correspond to the boundary values of the problem. 
%
% Parameters
% ==========
%    U - The2D matrix specifying the boundary conditions 
%        and which points to approximate the solution.
%
% Return Values
% =============
%    U_out - The matrix Uout is the 2D or 3D matrix which contains the 
%            approximation o the solution to Laplace’s equation. 

function [U_out] = laplace2d( U )
    [n_x, n_y] = size( U );
    
    % Error and Warning Checking
    % ==========================
    %
    % Test validaty of argument being passed in 
    % loop through all boundary values and check to ensure its not -Inf

    % ensure matrix dimentions are valid
    if n_x < 2 || n_y < 2
         throw( MException( 'MATLAB:invalid_argument', ...
            'the matrix passed in is not 2D' ) );
    end
    %ensure all boundry points are not -Inf
    for n = 1:n_x
        if U(n,1) == -Inf || U(n,end) == -Inf
            throw( MException( 'MATLAB:invalid_argument', ...
            'the matrix has a boundry value of -Inf' ) );
        end
    end
    for j = 1:n_y
        if U(1,j) == -Inf || U(end,j) == -Inf
             throw( MException( 'MATLAB:invalid_argument', ...
            'the matrix has a boundry value of -Inf' ) );
        end
    end
    
    % Initialization
    % ==============
    %
    %   Make a copy of the input matrix
    
    U_out = U;

    % Mapping the unknown points to a unique number from 1 to m
    % =========================================================
    %
    %   loop through all interior points
    %   for each point that is -Inf map coordinate to some unique number 
    %   between 1 and m

    % u_to_w specifies the ith unknown point which corresponds to the 
    % ith column in w_to_u
    u_to_w = zeros( n_x, n_y );
    % w_to_u specifies the coodinates of the ith unknown value in ith 
    % column
    w_to_u = zeros( 2, n_x * n_y );
    m = 0;

    for ix = 1:n_x
        for iy = 1:n_y
            if U(ix, iy) == -Inf
                m = m + 1;
                u_to_w(ix, iy) = m;
                w_to_u(:, m) = [ix, iy]';
            end
        end
    end

    
    % Creating and solving a system of linear equations
    % =================================================
    %
    %   for each unknown point, update the sparse matrix and vector appropriately. 
    %   Solve the system to get appromimations to the solution at each unknown point

    % Create the sparse system of linear equations
    M = spalloc( m, m, 5*m );
    b = zeros( m, 1 );

    for k = 1:m
        % Get the coordinates of the kth point
        % For each of the 4 adjacent points, determine if
        % the point is an insluated boundary point, a Dirichlet
        % boundary point or an unknown value and modify M as appropriate.
        c = w_to_u(:,k);
        coord_1 = c + [-1  0]';
        p_1 = U(coord_1(1),coord_1(2));
        coord_2 = c + [ 1  0]';
        p_2 = U(coord_2(1),coord_2(2));
        coord_3 = c + [ 0 -1]';
        p_3 = U(coord_3(1),coord_3(2));
        coord_4 = c + [ 0  1]';
        p_4 = U(coord_4(1),coord_4(2));
        if p_1 ~= -Inf && ~isnan(p_1)%Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p_1;
        elseif p_1 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(coord_1(1),coord_1(2));
            M(k,j) = M(k,j) + 1;
        end
        if p_2 ~= -Inf && ~isnan(p_2) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p_2;
        elseif p_2 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(coord_2(1),coord_2(2));
            M(k,j) = M(k,j) + 1;
        end    
        if p_3 ~= -Inf && ~isnan(p_3) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p_3;
        elseif p_3 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(coord_3(1),coord_3(2));
            M(k,j) = M(k,j) + 1;
        end
        if p_4 ~= -Inf && ~isnan(p_4) %Dirichlet condition
            M(k,k) = M(k,k) -1;
            b(k) = b(k) - p_4;
        elseif p_4 == -Inf %another unknown
            M(k,k) = M(k,k) -1;
            j = u_to_w(coord_4(1),coord_4(2));
            M(k,j) = M(k,j) + 1;
        end
    end

    w = M \ b;

    % Substituting the values back into the matrix U_out
    % ===================================================
    %
    %   mapping back from w space to u space.

    for k = 1:m
        % Copy the value from w into U_out
        point = w_to_u(:,k);
        value = w(k);
        U_out(point(1),point(2)) = value;
    end
end
