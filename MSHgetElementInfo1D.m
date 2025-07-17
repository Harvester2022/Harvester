function [long , rc , b ] = MSHgetElementInfo1D( obj )
%MSHgetElementInfo1D( obj ): Returns information of 1D (curves) domains in the mesh. 
%INPUT
	%	- obj: femSPACE2D object in which the mesh is embedded. 
	%OUTPUT: 
	% 	- long: nelement x 1 columne vector. Length of each line element. 
	%	- rc: nelement x 3 array. XYZ coordinates of mass center of each line element
    %	- b: nelement x 1 columne vector. Gradient of each element.

   
    %Load the mesh
    elementos = obj.msh.LINES;

    %Node coordinates
        x = obj.msh.POS1(:,1);
        y = obj.msh.POS1(:,2);
        z = obj.msh.POS1(:,3);

    %Position of 1 and 2 nodes in each elements
        x1 = x( elementos(:,1) );
        x2 = x( elementos(:,2) );

        y1 = y( elementos(:,1) );
        y2 = y( elementos(:,2) );

        z1 = z( elementos(:,1) );
        z2 = z( elementos(:,2) );

    %Relaticve positions
        %12 
        x12 = x2-x1;
        y12 = y2-y1;
        z12 = z2-z1;
        
    %Mass center coordinates
        rc( : , 1 ) = x1 + ( x12 )/2;
        rc( : , 2 ) = y1 + ( y12 )/2;
        rc( : , 3 ) = z1 + ( z12 )/2;       

        
    %Lenght of each elemen    
        long = sqrt( x12.^2 + y12.^2 + z12.^2 );
        
    %The gradient
        b1 = -long.^(-1);
        b2 = +long.^(-1);
        b = [ b1 , b2 ]; 