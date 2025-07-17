function F = CREATE_Js_Mag2D( obj , NIs , uM_P0 )
%ENSEMBLE_Js_Mag2D: Creates a load vector corresponding to the magnetic currents inside a NI surface domain with a given magnetization in XY plane. 
%INPUT
%	- NI: integer vector. Surface domain identifier. If 0, computed in all domains
%   - uM_P0: bidimensional P0 vector field with the XY plane magnetization. 

%OUTPUT: 
% 	- F. corresponding nbTRIANGLE x 1 element vector. 

        tic;

        uMx_P0 = obj.FIELD_component( uM_P0 , 'x' );
        uMy_P0 = obj.FIELD_component( uM_P0 , 'y' );

        [ielementos] = MSHgetElementsNodes( obj , 2 , NIs );

        TRI = obj.msh.TRIANGLES( ielementos , : );

        %Coordenates nodes in each TRiangle
            x1 = obj.msh.POS(TRI(:,1),1);
            x2 = obj.msh.POS(TRI(:,2),1);
            x3 = obj.msh.POS(TRI(:,3),1);
            y1 = obj.msh.POS(TRI(:,1),2);
            y2 = obj.msh.POS(TRI(:,2),2);
            y3 = obj.msh.POS(TRI(:,3),2);

        %Relativ position in each triangle
            Dx( : , 1 ) = x2-x1;
            Dx( : , 2 ) = x3-x2;
            Dx( : , 3 ) = x1-x3;
            Dy( : , 1 ) = y2-y1;
            Dy( : , 2 ) = y3-y2;
            Dy( : , 3 ) = y1-y3;


        %lenght of each edge in a triangle
            Long = sqrt( Dx.^2 + Dy.^2 );

        %Normal vectors
            nx = +Dy./Long;
            ny = -Dx./Long;           

            F = sparse( obj.msh.nbNod , 1 );

            nodos1 = TRI(:,[1,2,3]);
            nodos2 = TRI(:,[2,3,1]);

        for ie = 1:3                
            %magnetization currents in each triangle
            JmagP0 = uMx_P0( ielementos ).*ny( : , ie ) - uMy_P0( ielementos ).*nx( : , ie );

            %lenght of each edge in a triangle
            L = Long( : , ie );

            %Each edge has a contribution
            F = F + sparse( nodos1(:,ie) , 1 , JmagP0.*L/2 , obj.msh.nbNod , 1 );
            F = F + sparse( nodos2(:,ie) , 1 , JmagP0.*L/2 , obj.msh.nbNod , 1 );                
        end
        
        CPU = toc;
        if obj.verbosity == 1
            fprintf( 'Vector Jmag ensembled in %f s.\n' , CPU );      
        end
end