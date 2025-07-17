function [ elementos , nodos ] = MSHgetElementsNodes( obj , dim , NIDvec )
%MSHgetElementsNodes: get the element list and node list of domains in
%NIDvec integer vector
%   INPUT: 
%       - obj. femSPACE3D class object
%       - dim. Dimension of the domain. 0 or ‘Point’ for a point (0D). 1 or
%       ‘Line’ for a surface (1D). 2 or ‘SurfaceIntegral’ for a surface(2D).
%       3 or `Volume´ for a volume
%       - NIDvec. a vector with the domain ID numbers (integers)


    elementos0 = [];
    nodos0     = [];
  
    for iNID = 1:length( NIDvec )
        
        NID = NIDvec( iNID );

        switch dim
            
            case {0,'Point'}
                %caso punto
                nT = size( obj.msh.PNT , 1 );   

                %los nodos y elementos en el dominio
                listaElementos = 1:1:nT;

                boleana1 = obj.msh.LINES(:,2)==NID;
                elementos = listaElementos( boleana1 )';
                nodos     = unique( obj.msh.PNT(boleana1,1:1) );
                
                if NID==0                    
                    elementos = listaElementos';
                    nodos     = unique( obj.msh.PNT(:,1:1) );                    
                end
                
            case {1,'Line'}
                %caso curva
                nT = size( obj.msh.LINES , 1 );   

                %los nodos y elementos en el dominio
                listaElementos = 1:1:nT;

                boleana1 = obj.msh.LINES(:,3)==NID;
                elementos = listaElementos( boleana1 )';
                nodos     = unique( obj.msh.LINES(boleana1,1:2) );
                
                if NID==0                    
                    elementos = listaElementos';
                    nodos     = unique( obj.msh.LINES(:,1:2) );                    
                end
                
            case {2,'Surface'}
                %caso superficie
                
                nT = size( obj.msh.TRIANGLES , 1 );

                %los nodos y elementos en el dominio
                listaElementos = 1:1:nT;

                boleana1 = obj.msh.TRIANGLES(:,4)==NID;
                elementos = listaElementos( boleana1 )';
                nodos     = unique( obj.msh.TRIANGLES(boleana1,1:3) );                
                if NID==0                    
                    elementos = listaElementos';
                    nodos     = unique( obj.msh.TRIANGLES(:,1:3) );                    
                end
                
                
            case {3,'Volume'}
                %caso volumenes
                
                nT = size( obj.msh.TETS , 1 );

                %los nodos y elementos en el dominio
                listaElementos = 1:1:nT;

                boleana1 = obj.msh.TETS(:,5)==NID;
                elementos = listaElementos( boleana1 )';
                nodos     = unique( obj.msh.TETS(boleana1,1:4) );                
                if NID==0                    
                    elementos = listaElementos';
                    nodos     = unique( obj.msh.TETS(:,1:4) );                    
                end
                
        end
        
        elementos0 = [elementos0;elementos];
        nodos0 = [nodos0;nodos];
        
    end
    
    %Los nodos tienen que ser unicos
        nodos = unique( nodos0 );
        elementos = elementos0;


        

