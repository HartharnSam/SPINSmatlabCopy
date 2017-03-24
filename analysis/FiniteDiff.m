function [ varargout ] = FiniteDiff(x, OrdD, n, varargin)
% FiniteDiff : Create a finite difference matrix of arbitrary order for an
% arbitrary grid. Note that the algorithm will produce a differentiation
% with accuracy of at least order n, but it might be higher.
%
% x       -> Grid along which to compute derivatives.
%         ->    If x is a vector, then a 1D derivative is computed. If x is
%               cell array of vectors, then the N-D derivatives are
%               returned, each of which is created using kron. In this
%               case, sparsity is forced.
% OrdD    -> Order of differentiation (0 returns identity matrix)
% n       -> Order of accuracy
% sp      -> Is sparse? (bool) (default: false)
% uniform -> Is unifrom grid? (bool) (default: true)
%
% There is a known bug, that odd order (OrdD) greater than 1 returns garbage.
% Just use even orders.
%
% courtesy Benjamin Storer, 2015

if length(varargin) == 0
    sp = false;
    uniform = true;
elseif length(varargin) == 1
    sp = varargin{1};
    uniform = true;
elseif length(varargin) == 2
    sp = varargin{1};
    uniform = varargin{2};
elseif length(varargin) > 2
    error('Too many optional arguments given. Two are allowed: Sparsity and Uniformity')
end

if iscell(x)
    % Create the differentiation matrix for each dimension's grid.
    numDimensions = length(x);
    Dx = cell(1,numDimensions);
    for ii = 1:numDimensions
        Dx{ii} = MakeMatr(x{ii}, n, true, uniform, OrdD);
    end
    % Now we do all of the krons in order to make them apply to a
    % vectorized grid
    for ii = 1:numDimensions % Loop through each dimension
        I = speye(length(x{ii}));
        for jj = 1:numDimensions % Loop through each derivative
            
            % Note: Since meshgrid puts x in the second dimension and y in
            % the first dimension, we deal with those two cases seperately.
            if (ii == 1) && (jj == 2)
                Dx{jj} = kron(I,Dx{jj});
            elseif (ii == 2) && (jj == 1)
                Dx{jj} = kron(Dx{jj},I);
            else
                if jj < ii
                   Dx{jj} = kron(I,Dx{jj});
                elseif jj > ii
                   Dx{jj} = kron(Dx{jj},I);
                end
            end
        end
    end
else
    Dx = MakeMatr(x, n, sp, uniform, OrdD);
    numDimensions = 1;
end


if numDimensions > 1
    for ii = 1:nargout
        varargout{ii} = Dx{ii}; %#ok<AGROW>
    end
else
    varargout{1} = Dx;
end

end

function Dx = MakeMatr(x, n, sp, uniform, OrdD)

if length(x) == 3
    % Using a length of 3 is shorthand. x = [a,b,c] is interpreted as
    % x = linspace(a,b,c). The advantage is that we don't actually need
    % to generate x, so we can save on memory when we want large grids.
    Nx = x(3);
else
    Nx = length(x);
end

% Adjustment for higher order derivatives
n = n + OrdD - 1;

if OrdD == 0
    Dx = eye(length(x));
    return;
end

warning('off', 'MATLAB:nearlySingularMatrix')

if sp
    Dx = sparse(Nx, Nx);
else
    Dx = zeros(Nx, Nx);
end

if uniform
    if length(x) == 3
        % Using a length of 3 is shorthand. x = [a,b,c] is interpreted as
        % x = linspace(a,b,c). The advantage is that we don't actually need
        % to generate x, so we can save on memory when we want large grids.
        dx = (x(2) - x(1))/(x(3)-1);
    else
        dx = x(2) - x(1);
    end
    % Deal with boundary issues
    for i = 1:ceil(n/2)
        A = zeros(n+1,n+1);
        for j = 1:n+1
            A(:,j) = (((j-i)*dx).^(0:n))./factorial(0:n);
        end
%         display(A)
        b = zeros(n+1,1);
        b(1+OrdD) = 1;
%         display(b)
        coeff = A\b;
        coeff = coeff';
        Dx(i, 1:n+1) = coeff;
    end
    for i = Nx-ceil(n/2)+1:Nx
        A = zeros(n+1,n+1);
        for j = Nx-n:Nx
            A(:,j-Nx+n+1) = (((j-i)*dx).^(0:n))./factorial(0:n);
        end
        b = zeros(n+1,1);
        b(1+OrdD) = 1;
        coeff = A\b;
        coeff = coeff';
        Dx(i, Nx-n:Nx) = coeff;
    end
    % Now do the internals.
    A = zeros(n+1,n+1);
    if mod(n,2) == 0 % If even...
        for j = -n/2:n/2
            A(:,j+n/2+1) = ((j*dx).^(0:n))./factorial(0:n);
        end
        b = zeros(n+1,1);
        b(1+OrdD) = 1;
        coeff = A\b;
        coeff = coeff';
        coeff = repmat(coeff, Nx-n, 1);
        Dx(n/2+1:Nx-n/2,:) = spdiags(coeff, 0:n, Nx-n, Nx);
    elseif mod(n,2) == 1 % If odd...
        for j = -floor(n/2):ceil(n/2)
            A(:,j+floor(n/2)+1) = ((j*dx).^(0:n))./factorial(0:n);
        end
        b = zeros(n+1,1);
        b(1+OrdD) = 1;
        coeff = A\b;
        coeff = coeff';
        coeff = repmat(coeff, Nx-n, 1);
        if OrdD == 1
            Dx(ceil(n/2)+1:Nx-ceil(n/2),:) = spdiags(coeff, 0:n, Nx-n-1, Nx);
        elseif OrdD > 1
            Dx(ceil(n/2):Nx-ceil(n/2),:) = spdiags(coeff, 0:n, Nx-n, Nx);
        end
    end
else
    for i = 1:Nx
        if i <= ceil(n/2)
            % Deal with boundary issues
            A = zeros(n+1,n+1);
            for j = 1:n+1
                dx = x(j)-x(i);
                A(:,j) = (dx.^(0:n))./factorial(0:n);
            end
            b = zeros(n+1,1);
            b(1+OrdD) = 1;
            coeff = A\b;
            coeff = coeff';
            Dx(i, 1:n+1) = coeff;
        elseif i > Nx - n
            % Deal with boundary issues
            A = zeros(n+1,n+1);
            for j = Nx-n:Nx
                dx = x(j)-x(i);
                A(:,j-Nx+n+1) = (dx.^(0:n))./factorial(0:n);
            end
            b = zeros(n+1,1);
            b(1+OrdD) = 1;
            coeff = A\b;
            coeff = coeff';
            Dx(i, Nx-n:Nx) = coeff;
        elseif i <= (Nx-ceil(n/2))
            % Deal with the internal pieces
            % If n is even, then just use a centred scheme
            if mod(n,2) == 0
                A = zeros(n+1,n+1);
                for j = -n/2:n/2
                    dx = x(i+j) - x(i);
                    A(:,j+n/2+1) = (dx.^(0:n))./factorial(0:n);
                end
                b = zeros(n+1,1);
                b(1+OrdD) = 1;
                coeff = A\b;
                coeff = coeff';
                Dx(i, i-n/2:i+n/2) = coeff;

            % If n is odd, then bias to which side has the closest point.
            elseif mod(n,2) == 1
                if abs(x(i+ceil(n/2)) - x(i)) <= abs(x(i-ceil(n/2)) - x(i))
                    A = zeros(n+1,n+1);
                    for j = -floor(n/2):ceil(n/2)
                        dx = x(i+j) - x(i);
                        A(:,j+floor(n/2)+1) = (dx.^(0:n))./factorial(0:n);
                    end
                    b = zeros(n+1,1);
                    b(1+OrdD) = 1;
                    coeff = A\b;
                    coeff = coeff';
                    Dx(i, i-floor(n/2):i+ceil(n/2)) = coeff;
                else
                    A = zeros(n+1,n+1);
                    for j = -ceil(n/2):floor(n/2)
                        dx = x(i+j) - x(i);
                        A(:,j+ceil(n/2)+1) = (dx.^(0:n))./factorial(0:n);
                    end
                    b = zeros(n+1,1);
                    b(1+OrdD) = 1;
                    coeff = A\b;
                    coeff = coeff';
                    Dx(i, i-ceil(n/2):i+floor(n/2)) = coeff;
                end
            end
        end
    end
end

warning('on', 'MATLAB:nearlySingularMatrix')

end
