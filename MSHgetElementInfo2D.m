function [Area , rc , un , b , c ] = MSHgetElementInfo2D( obj )
%MSHgetElementInfo2D( obj ): Returns information of 2D (surface) domains in the mesh. 
%INPUT
	%	- obj: femSPACE2D object in which the mesh is embedded. 
	%OUTPUT: 
	% 	- Area: nelement x 1 column vector. Area of each triangle element. 
	%	- rc: nelement x 3 array. XYZ coordinates of mass center of each triangle. 
	%	- un: P0 vector field defined in each 2D surface. It corresponds to the unitary vector perpendicular to each triangle (surface).
  
   
    %Load the mesh
    TRI = obj.msh.TRIANGLES;
    
    % nodes 1 2 3 of each triangle corners
        nodos1 = TRI( : , 1 );
        nodos2 = TRI( : , 2 );
        nodos3 = TRI( : , 3 );

    %Positions
        x = obj.msh.POS1(:,1);
        y = obj.msh.POS1(:,2);
        z = obj.msh.POS1(:,3);
        
    % Position of each node in the triangle
        x1 = x( nodos1 );
        x2 = x( nodos2 );
        x3 = x( nodos3 );
        
        y1 = y( nodos1 );
        y2 = y( nodos2 );
        y3 = y( nodos3 );
        
        z1 = z( nodos1 );
        z2 = z( nodos2 );
        z3 = z( nodos3 );

    %Relative postion vectores in each triangle
        r1 = [ x1 , y1 , z1 ];
        r2 = [ x2 , y2 , z2 ];
        r3 = [ x3 , y3 , z3 ];

        r12 = r2-r1;
        r13 = r3-r1;        
        
    %Cross product for the perpendicular vector. 
        nvec = cross( r12 , r13 )./sqrt( sum( cross( r12 , r13 ).^2 , 2 ) );
        
    %MC positions
        rc = r1 + ( r12 + r13 )/2;
        
    %The perpendicular unitary vector field
        un = reshape( nvec , [] , 1 );
        

        
    %Areas of each triangle element     
        Area = ( sqrt( sum( cross( r12 , r13 ).^2 , 2 ) ) )/2;
        
    %The gradients
        b1 = ( y2-y3 )./Area/2;
        b2 = ( y3-y1 )./Area/2;
        b3 = ( y1-y2 )./Area/2;
        c1 = ( x3-x2 )./Area/2;
        c2 = ( x1-x3 )./Area/2;
        c3 = ( x2-x1 )./Area/2;
        
        b = [ b1, b2 , b3 ];
        c = [ c1, c2 , c3 ];
 