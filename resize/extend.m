function array3D = extend(array2D,P,dim)
%EXTEND Extends a 2D array in a third dimension. 
%   
%           array2D       : An MxN array
%           N (int)       : Number of pages in third dimension
%           dimension     : Sets P in the dimension specified by dim.
%
%           Options       : dim = 1: size(array3D) = PxMxN 
%                           dim = 2: size(array3D) = MxPxN
%                           dim = 3: size(array3D) = MxNxP


    sz = size(array2D);
    Nx = sz(1);
    Nz = sz(2);
    
    M = zeros(Nx,Nz,1);
    
    M(:,:,1) = array2D;
    M = repmat(M,[1 1 P]);
    %size(M)
    if dim == 1
        M = permute(M,[3 1 2]);
    elseif dim == 2
        M = permute(M,[1 3 2]);
    elseif dim == 3
        M = permute(M,[1 2 3]);
    else
        error('new_dim must be 1, 2 or 3')
    end
    %size(M)
    array3D = M;
   
end

