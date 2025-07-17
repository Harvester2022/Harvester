function S = MSHgetDomainSices( obj , dim , NIDvec )
%MSHgetElementsSices: get the sizes (lenght, area, volume) of the element of domains in
%NIDvec integer vector
%   INPUT: 
%       - obj. femSPACE3D class object
%       - dim. Dimension of the domain. 0 or ‘Point’ for a point (0D). 1 or
%       ‘Line’ for a surface (1D). 2 or ‘SurfaceIntegral’ for a surface(2D).
%       3 or `Volume´ for a volume
%       - NIDvec. a vector with the domain ID numbers (integers)


   elementos = MSHgetElementsNodes( obj , dim , NIDvec );
   
    switch dim           
                
        case {1,'Line'}
            %caso curva
            usize = obj.msh.D1.longitude( elementos );

        case {2,'Surface'}
            %caso superficie                
            usize = obj.msh.D2.Area( elementos );               

        case {3,'Volume'}
            %caso volumenes
            usize = obj.msh.D3.Volume( elementos );                  
    end
    
    S = sum( usize );