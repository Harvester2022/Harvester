classdef femSPACE2D <handle
    %femSPACE2D: This is a MatLab class for the implementation of 1st order 2D finite elemetns. For a provided mesh, it creates the corresponding fem matrix and vectors from the variation formulation. The mesh must be created from gmsh.
    properties
        %u field dimensions: 1 or 2. Default 1D (escalar). 2D (vectorial)
        udim = 1;
        %Scale. The geometry is resized by escala. Useful when changin lenght units. 
        scale = 1;
        %Problem name
        name = 'fem2D';
        %verbosity (1 = activate output messages, default)
        verbosity = 1;
        %struc to save mesh information
        msh;
    end

    methods
        
        function obj = femSPACE2D( gmshFile , udim , scale , ParameterList )
        %femSPACE2D: Reads the mesh file from gmsh and load it in femSPACE2D object for the finite element implementation. It applies a scale to reshide the mesh.
            % INPUT:
            %	- gmshFile: gmsh file where the mesh. If the file is .m (mesh to be read by matlab) it only reads an load it. If the file is .geo (gmsh input), it calls gmsh to perform the meshing and create the tempMSH.m.m file with the mesh. Afterwards it loads tempMSH.m file.
            %   - udim: dimensions of the unknown field. Default 1
            %	- scale: scale factor. Default 1, no scale. 


            if ~exist('scale','var')
                 % default value
                  scale = 1;
             end
             
             obj.scale = scale;
             
             if exist('ParameterList','var')
                 fileID = fopen('Parameters.geo','w');
                 NameList = fieldnames( ParameterList );
                 for i1 = 1:length( NameList )
                     Pname = NameList{i1};
                     Pvalue = getfield( ParameterList , Pname );
                     fprintf(fileID,'%s = %12e;\n',Pname,Pvalue);
                 end
                 fclose(fileID);
             end
             
             %u field dimensions
             obj.udim = udim;

             %load gmsh.m file
             %if it is a .geo field
             if gmshFile(end-3:end)=='.geo'

                 %mesh, create tempMSH.m file and load it
                 comado1 = sprintf( 'gmsh.exe -3 %s -format m -o ./tempMSH.m' , gmshFile );
                 system( comado1 );
                 tempMSH;
                 system( 'del tempMSH.m' );

             else
                %Si no, se carga
                eval( gmshFile );
             end

            %rescale
            POS1 = msh.POS*scale;
            msh.POS  = POS1;
            msh.POS1 = POS1;

            %save the mesh information
            obj.msh = msh;
            
            %call prepareMesh to finish the calculations
            obj.prepareMesh( );
            
       end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
        function prepareMesh( obj )
        %prepareMesh: performs geometric calculations necessary to ensemble the matrix. Run every time a new mesh is loaded. 
            %Elementos de 1D. Si existen
            if isfield( obj.msh , 'LINES' )==1
                [longitude , rc1 , b ] = MSHgetElementInfo1D( obj );
                obj.msh.D1.longitude = longitude;
                obj.msh.D1.rc = rc1;
                %1D gradient
                obj.msh.D1.Grad.b = b;
                %Line number
                obj.msh.nbLINE = size( obj.msh.LINES  , 1 );          
            end

            %Elementos de 2D. Si eisten
            if isfield( obj.msh , 'TRIANGLES' )==1
                [~ , ~ , nvec , ~ , ~ ] = MSHgetElementInfo2D( obj );
                
                %if the z component of nvec is negative, permute triangle
                %node numbers in order to be positive, e.g, change 1 with 2. 
                nvec = reshape( nvec , [] , 3 );
                obj.msh.TRIANGLES( nvec(:,3)<0 , [1,2] ) = obj.msh.TRIANGLES( nvec(:,3)<0 , [2,1] );
                %recalculate after the permutation                    
                [Area , rc2 , nvec , b , c ] = MSHgetElementInfo2D( obj );
                obj.msh.D2.Area = Area;
                obj.msh.D2.rc = rc2;
                %triangle normal vector
                obj.msh.D2.nvec = nvec;
                %triangle number
                obj.msh.nbTRI = size( obj.msh.TRIANGLES  , 1 );
                %2D gradient
                obj.msh.D2.Grad.b = b;
                obj.msh.D2.Grad.c = c;                
            end
        end



        function M = createMatrix( obj , domainDimensions , NI , type , fP0 )
        %createMatrix: creates a fem matrix for the field calculations. The matrix is equivalent to the weak formulation described in type and integrated acros NI domain of DomainDimensions. 
        %INPUT
        %	- DomainDimensions: dimensions of the domain where the integral
        %	is performed. 1 or ‘LineIntegral’ for a curve (1D). 2 or
        %	‘SurfaceIntegral’ for a surface(2D). 3 or 'VolumeIntegral' for
        %	a volume (3D)
        %   - NI: integer vector. Identifier of the domains where the the integral is performed. If 0 is provided, the integral will cover all the domains with DomainDimensions dimnensionality. 
        %	- type: String. Body of the integral in variational form. See documentation. 
        %	- fP0: scalar field defined in each element (P0). If a number is provided, f will be constant in the NI domanis. 
        %OUTPUT: 
        %	- M: sparse matrix. fem matrix for field calculations


            %load node number
            nNodos = obj.msh.nbNod;
            %initialise the matrix
            M = sparse( nNodos , nNodos );

            %si no se da fP0, es uno
            if ~exist('fP0','var')
                fP0=1;
            end
            
            switch domainDimensions
                


                %in 2 dimension.........................................
                case {2,'SurfaceIntegral'}

                    [ ielementos , ~] = MSHgetElementsNodes( obj , 2 , NI );
                    Areas = obj.msh.D2.Area( ielementos );
                    b = obj.msh.D2.Grad.b( ielementos , : );
                    c = obj.msh.D2.Grad.c( ielementos , : );
                    
                    switch type
                        case {'dv/dx*du/dx' , '1,1'}
                            g = @(i,j)b(:,i).*b(:,j);
                        case {'dv/dx*du/dy' , '1,2'}
                            g = @(i,j)b(:,i).*c(:,j);
                        case {'dv/dy*du/dx' , '2,1'}
                            g = @(i,j)c(:,i).*b(:,j);
                        case {'dv/dy*du/dy' , '2,2'}
                            g = @(i,j)c(:,i).*c(:,j);
                        case {'v*u' , '0,0'}
                            g = @(i,j) ( (i==j) + 1 )/12;
                        case {'v*du/dx' , 'D1'}
                            g = @(i,j)b(:,j)/3;
                        case {'v*du/dy' , 'D2'}
                            g = @(i,j)c(:,j)/3;                            
                        otherwise
                        fprintf('createMatrix. SurfaceIntegral. ERROR! Incorrect variational form! \n')
                    end

                    %if fP0 is a P0 field
                    if length( fP0 )>1
                        fP0 = fP0( ielementos );
                    end

                    for i=1:3
                        for j=1:3
                            %for each node in the triangle
                            indxtri_i = obj.msh.TRIANGLES(ielementos,i);
                            indxtri_j = obj.msh.TRIANGLES(ielementos,j);
                            M = M + sparse( indxtri_i , indxtri_j , g(i,j).*fP0.*Areas , nNodos , nNodos );
                        end
                    end


                %in 1 dimension.........................................
                case {1,'LineIntegral'}

                    [ ielementos , ~ ] = MSHgetElementsNodes( obj , 1 , NI );
                    L = obj.msh.D1.longitude( ielementos );

                    switch type
                        case 'v*u'
                            g = @(i,j) ( (i==j) + 1 )/6;
                        otherwise
                            fprintf('createMatrix. LineIntegral. ERROR! Incorrect variational form! \n')
                    end

                    %if fP0 is a P0 field
                    if length( fP0 )>1
                        fP0 = fP0( ielementos );
                    end

                    for i=1:2
                        for j=1:2
                            %for each node in the LINE
                            indxtri_i = obj.msh.LINES(ielementos,i);
                            indxtri_j = obj.msh.LINES(ielementos,j);
                            M = M + sparse( indxtri_i , indxtri_j , g(i,j).*fP0.*L , nNodos , nNodos );
                        end
                    end
                otherwise 
                    fprintf('createMatrix. ERROR! Incorrect domainDimensions! \n')
            end
        end


        function F = createVector( obj , domainDimensions , NI , fP0 )
        %createVector: creates a fem vector for the field calculations. The vector is equivalent to the weak formulation integrated across NI domain of DomainDimensions.
        %INPUT
        %	- DomainDimensions: dimensions of the domain where the integral is performed. 1 or ‘LineIntegral’ for a curve (1D). 2 or ‘SurfaceIntegral’ for a surface(2D)
        %   - NI: integer vector. Identifier of the domains where the the integral is performed. If 0 is provided, the integral will cover all the domains with DomainDimensions dimnensionality. 
        %	- fP0: scalar field defined in each element (P0). If a number is provided, f will be constant in the NI domanis. 

        %OUTPUT: 
        %	- F: sparse vector. fem vector for field calculations


            %load node number
            nNodos = obj.msh.nbNod;
            %initialize the vector
            F = sparse( nNodos , 1 );
            if ~exist('fP0','var')
                fP0=1;
            end

            switch domainDimensions

                case {2,'SurfaceIntegral'}
                    [ ielementos , ~] = MSHgetElementsNodes( obj , 2 , NI );
                    Area = obj.msh.D2.Area( ielementos );
                    %if fP0 is a P0 field
                    if length( fP0 )>1
                        fP0 = fP0( ielementos );
                    end
                    
                    for i=1:3
                        %for each node in the triangle
                        indxtri_i = obj.msh.TRIANGLES(ielementos,i);
                        F = F + sparse( indxtri_i , 1 , fP0.*Area/3 , nNodos , 1 );
                    end

                case {1,'LineIntegral'}
                    [ ielementos , ~] = MSHgetElementsNodes( obj , 1 , NI );
                    L = obj.msh.D1.longitude( ielementos );
                    %if fP0 is a P0 field
                    if length( fP0 )>1
                        fP0 = fP0( ielementos );
                    end
                    for i=1:2
                        %for each node in the LINE
                        indxtri_i = obj.msh.LINES(ielementos,i);
                        F = F + sparse( indxtri_i , 1 , fP0.*L/2 , nNodos , 1 );
                    end
                otherwise
                    fprintf('createVector. ERROR! Incorrect domainDimensions! \n')
            end
        end
        
        
        function M = create2DDerivativeP0( obj , NI , type , fP0 )
        %create2DDerivativeP0: creates a Dij matrix in such a way that f*du/dx is given by dudx = D*u. The resulting dudx derivative is a P0 field whereas u is a P1 field. 
        %INPUT
        %   - NI: integer vector. Identifier of the domains where the the derivative is performed. 
        %	- type.: ‘dx’ or 1, derivative over x. If ‘dy’ or 2, derivative over y. 
        %	- fP0: P0 field.  The values of f. If a single number is provided, f is assumed to be constant. 

        %OUTPUT: 
        %	- Dij: sparse nbTRI x nbNodes matrix. 

                %load
                nnodos = obj.msh.nbNod;
                nelementos = obj.msh.nbTRI;        
                [ ielementos , ~] = MSHgetElementsNodes( obj , 'Surface' , NI );
                
                if ~exist('fP0','var')
                    fP0=1;
                end
                
                if length( fP0 )>1
                    fP0 = fP0( ielementos );
                end
                
                M = sparse( nelementos , nnodos );
            
                switch type            
                    case {'dx',1}               
                        BC = obj.msh.D2.Grad.b;
                    case {'dy',2}                
                        BC = obj.msh.D2.Grad.c; 
                    otherwise
                        fprintf('Create2DDerivativeP0. ERROR! Incorrect variational form!')
                end    
                                
                for i1=1:3
                        %for each nodes in element
                        indxtri_i = obj.msh.TRIANGLES( ielementos , i1 );
                        M = M + sparse( ielementos , indxtri_i , BC( ielementos , i1 ).*fP0 , nelementos , nnodos );    
                end           
        end
        
        
        
        
         
        function A = ENSEMBLE_Laplace2D( obj , NI , fP0x , fP0y )
            %ENSEMBLE_Laplace2D: Ensembles the Stiffness matrix for surface domains identified by NI. The stiffness matrix correspond to the simple Laplace operator applied to the scalar u field. 
            %   Dx(fx*Dx(u)) + Dy(fy*Dy(u))
            %INPUT
            %	- obj: femSPACE2D object in which the mesh is embedded. 
            %	- NI: integer vector. Domain identifier. If 0, computed in all domains. 
            %   - fP0x: P0 field with constant value at each element in NI. Default, unity for all elements            
            %   - fP0y: P0 field with constant value at each element in NI. Default, unity for all elements
            %OUTPUT: 
                % 	- LP: the matrix.
            tic; 
                if ~exist('fP0x','var')
                    fP0x = 1;
                end

                if ~exist('fP0y','var')
                    fP0y = fP0x;
                end
            A = - obj.createMatrix( 'SurfaceIntegral', NI , 'dv/dx*du/dx' , fP0x ) - obj.createMatrix( 'SurfaceIntegral', NI , 'dv/dy*du/dy' , fP0y );
            CPU = toc;            
            if obj.verbosity == 1
                fprintf( 'Laplace 2D matrix ensembled in %f s.\n' , CPU );
            end
        end
        
        
        
        
        
        
        
        function A = ENSEMBLE_curlcurl2D( obj , NI , fP0x , fP0y )
            %ENSEMBLE_curlcurl2D: Ensembles the Stiffness matrix for surface domains identified by NI. The stiffness matrix correspond to corresponding to doble curl operator aplied accross plane XY to rot(f.rot(0,0,uz)) vector field. 
            %   Dx(fx*Dx(u)) + Dy(fy*Dy(u))
            %INPUT
            %	- obj: femSPACE2D object in which the mesh is embedded. 
            %	- NI: integer vector. Domain identifier. If 0, computed in all domains. 
            %   - fP0x: P0 field with constant value at each element in NI. Default, unity for all elements            
            %   - fP0y: P0 field with constant value at each element in NI. Default, unity for all elements
            %OUTPUT: 
                % 	- LP: the matrix.
            tic; 
                if ~exist('fP0x','var')
                    fP0x = 1;
                end

                if ~exist('fP0y','var')
                    fP0y = fP0x;
                end
            A = +obj.createMatrix( 'SurfaceIntegral', NI , 'dv/dx*du/dx' , fP0y ) +obj.createMatrix( 'SurfaceIntegral', NI , 'dv/dy*du/dy' , fP0x );
            CPU = toc;            
            if obj.verbosity == 1
                fprintf( 'rotot 2D matrix ensembled in %f s.\n' , CPU );
            end
        end
        
        
        
        
        
        
        
        function Mmass = ENSEMBLE_Mass2D( obj , dim , NI , fP0 )
            %ENSEMBLE_Mass2D: Ensembles the mass matrix for surface domains identified by NI. 
            %INPUT
            %	- obj: femSPACE2D object in which the mesh is embedded. 
            %   - dim: deimension of NI domains
            %	- NI: integer vector. Domain identifiers. If 0, computed in all domains. 
            %   - fP0: P0 field with constant value at each element in NI. Default, unity for all elements

            %OUTPUT: 
                % 	- Mmass: the mass matrix. 

               tic;   

               %u field dimensions
                    uDim = obj.udim;         

               if ~exist('fP0','var')
                    fP0 = 1;
               end 

               Mmass = kron( eye( uDim ) , obj.createMatrix( dim , NI , 'v*u' , fP0 ) );

               CPU = toc;
               if obj.verbosity == 1
                fprintf( 'Mass matrix ensembled in %f s.\n' , CPU );
               end
        end
            
        function MS = ENSEMBLE_Stiffnes2D( obj , NI , Cijkl , ubeta )
            %ENSEMBLE_Stiffness2D: Ensembles the general stiffness matrix for surface domains identified by NI. 
            %INPUT
            %	- NI: integer vector. Domain identifier. If 0, computed in all domains. 
            %   - Cijkl: Cell(2x2x2x2). Each element C{i,j,k,l} is a P0 field with constant value (Cijkl) at each triangle in NI. 
            %   - ubeta: Rayleight damping coefficeint If not provided, 1, no effect. 
            %OUTPUT: 
            % 	- MS: the stiffness matrix. 
                tic;
                
               if ~exist('ubeta','var')
                    ubeta = 1;
               end 

                %NI element dimensions, surface: 
                dim = 2;
                nn = obj.msh.nbNod;
                [ni,nj,nk,nl] = size( Cijkl );
                MS = sparse( nn*ni , nn*ni );
                for i = 1:ni
                    %row indexs
                    ii = (1:nn) + (i-1)*nn;   
                    for j = 1:nj
                        %column indexs
                        jj = (1:nn) + (j-1)*nn;
                        for k=1:nk
                            for l=1:nl
                                if size( Cijkl{i,j,k,l} )== 0
                                    Cijkl{i,j,k,l}=0;
                                else
                                    tipoVariacional = sprintf( '%d,%d' , k , l );
                                    MS(ii,jj) = MS(ii,jj) + obj.createMatrix( dim , NI , tipoVariacional , Cijkl{i,j,k,l}.*ubeta );
                                end
                            end
                        end
                    end
                end
                CPU = toc;
                if obj.verbosity == 1
                    fprintf( '2D stiffness matrix ensambled in %d en %f s.\n' , NI(1) , CPU );  
                end
                
                %negative
                MS = -MS;
        end      
        
        
        
        
        
        function ML = ENSEMBLE_Load( obj , dim , NI , fP0 )
            %ENSEMBLE_Load: Ensembles the load vector for 1 (curve), 2 (surface) or 3 (volume) dimensional domains identified by NI  
            %INPUT
            %	- obj: femSPACE2D object in which the mesh is embedded. 
            %	- dim. dimensions of NI domain. If it is a curve (1D), 1 or ‘LineIntegral’. If it is a surface, 2 or ‘SurfaceIntegral’. If it is a surface 3 or ‘Volume’. 
            %	- NI: integer vector. Domain identifier. If 0, computed in all domains. 
            %   - fP0: P0 escalar or vector field with constant value at each element in NI. Default, unity for all elements
            %OUTPUT: 
                % 	- ML: Load vector. 
                tic;
                %field dimensions
                Udim = obj.udim;  
                %if no  fP0 field provided, 1. 
                if ~exist('fP0','var')
                    fP0 = ones(1,Udim);
                end   

                %fP0 field components nelements x udim
                fP0 = reshape( fP0 , [] , Udim );
                ML = [];    
                for i1=1:Udim    
                    ML = [ML;obj.createVector( dim , NI , fP0(:,i1) )];
                end
                CPU = toc;
                if obj.verbosity == 1
                    fprintf( 'Load 2D Vectors ensebled for domain %d  in %f s.\n' , NI(1) , CPU );
                end
            end

        
        
        
               
        
             
        
        
        %OCTAVE/MATLAB SOLVER********************************
        function u = SOLVER( obj , A , B )
        %SOLVER: Solves the u P1 field for the A*u=B linear equation system. 
        %INPUT
        %	- A: nbNod x nBnod array matrix. First, null column/rows are removed. 
        %           - B: nbNod x 1 array column vector. 
        %OUTPUT: 
            % 	- u: P1 field. 

                tic;
                if (obj.verbosity==1); fprintf('Starting Solver ... '  ); end
                
                %remove null submatrix
                inulos = any( A , 2 );
                a1 = A( inulos , inulos );
                b1 = B( inulos );

                u1 = a1\b1;
                
                u = zeros( size(B) );
                u( inulos ) = u1;
                
                u = full( u );

                CPU = toc;

                if (obj.verbosity==1); fprintf('Solver finished in %.3f s.\n' , CPU  );end

        end      
        
        
        
        
        
        function [v,u] = SOLVEmodes( obj , A , M , ID )
        %SOLVERmodes: Finds the autovalues and autovectors (modes) for the A*u=M*u linear equation system. First, null column/rows are removed. If ID is provided, homogeneous ditched (u=0) conditions will be applied at that curve domain. 
        %INPUT
        %	- A: nbNod x nBnod array matrix. Usually stiffness matrix. 
        %           - M: nbNod x nBnod array matrix. Usually mass matrix.
        %	- ID: (optional) integer vector with the ID numbers where homogeneous ditched (u=0) conditions will be applied. 
        %OUTPUT: 
            % 	- v: array with the first smaller autovalues.
            
                tic;                
                nvalues = 14;     
                udims = obj.udim;     
                %homogeneous Dirichlet
                if exist( 'ID' , 'var' )
                    A = A + kron( eye( udims ) , obj.createMatrix( 'LineIntegral' , ID , 'v*u' , 1e33 ) );    
                end               
                %remove null submatrix
                inulos = any( A , 2 );
                A1 = A( inulos , inulos );
                M1 = M( inulos , inulos );
                
                [u1,v1] = eigs( A1 , M1 , nvalues , 0 );
                v = diag( v1 );
                
                u = zeros( size(A,1) , nvalues );
                u( inulos , : ) = u1;
                
                CPU = toc;

                if (obj.verbosity==1);fprintf('Eigenvalues found in %.3f s.\n' , CPU  );end;

        end  

        
        
        
        
              
        function [A,B] = SOLVERapplyDirichlet( obj , A , B , dimE , NIB , valor , inDims )
        %SOLVERapplyDirichlet: Applies the Dirichlet conditions to boundary domains in ID.  The problem is A*u=B.
        %INPUT
        %	- A: nbNod x nBnod array matrix. 
        %           - B: nbNod x 1 array vector. 
        %	- dimelE: dimensions of NI domain. If it is a curve (1D), 1 or ‘Line’. If it is a surface, 2 or ‘Surface’. If it is a volume, 3 or ‘Volume’
        %	- ID: integer vector with the ID numbers where ditched conditions will be applied. 
        %OUTPUT: 
        %	- A: nbNod x nBnod array matrix with the Dirichlet boundary conditions applied.
        %	- B  nbNod x 1 array. array vector with the Dirichlet boundary conditions applied.
        
                %dimensions of u should be:
                dimu = obj.udim;
        
                if ~exist('inDims','var')
                    inDims = 1:1:dimu;  
                end 
                
                %boundary domain number
                nB = length( NIB ); 

                %node number
                nbNod = obj.msh.nbNod;

                %For each boundary domain
                for iB = 1:nB
                     [ ~ , nodos ] = MSHgetElementsNodes( obj , dimE , NIB );
                    for iD = 1:1:length( inDims )
                        %node index
                        indiceNodos = nodos + nbNod*( inDims(iD) - 1 );
                        A( indiceNodos ,      :  ) = 0;
                        A( indiceNodos , indiceNodos ) = eye( length( nodos ) );
                        B( indiceNodos ) = valor(iD);
                    end
                end                
        end    
        
        
        
        
        
        
        
        
        
        
        
        %POST PROCES FUNCTIONS********************************

        function [SIGMA] = POST_2DSIGMA_P0( obj , u , NI , Cijkl )
        %POST_2DSIGMA_P0: Calculates the sigma tensor defined as: sigma_i,k=C_i,j,k,l du_j/dx_l. where u is a P1 field. The resulting tensor is a P0 field and it is computed over a ID surface domains. 
        %INPUT
        %           - u: vector field P1
        %	- NI: Integer vector, Surface domains where sigma is computed. If 0, is computed over all domains. 
        %	- C_i,j,k,l: 2x2x2x2 cell.  Each cell element is a scalar P0 field whit the corresponding Cijkl coeficiens. If a number is provided, Cijkl is assumed to be constant across the domain. 

        %OUTPUT: 
        %	- SIGMA: sparse matrix array. Each colum is a bidimiensional vector P0 field. The column k correspond to sigma_i,k. 


            %dimensiones del campo
            [ni,nj,nk,nl] = size( Cijkl );

            %numero de nodos en el mallado
            nNodos = obj.msh.nbNod;
            
            %numero de tetraedros en el mallado
            nTRI = obj.msh.nbTRI;
            
            SIGMA = zeros( nTRI*ni , nk );

            for k = 1:nk
                %reinicial lor matrices
                Mgrad = sparse( nTRI*ni , nNodos*nj );
                for i = 1:ni
                    for j=1:nj
                        %indices
                        ii = (1:nTRI) + (i-1)*nTRI;
                        jj = (1:nNodos) + (j-1)*nNodos;
                        for l=1:nl
                            %Si no es elemento nulo calcular
                            if size( Cijkl{i,j,k,l} ) > 0
                                Mgrad(ii,jj) = Mgrad(ii,jj) + obj.create2DDerivativeP0( NI , l , Cijkl{i,j,k,l} );
                            end
                        end

                    end
                end
                SIGMA(:,k) = Mgrad*u;
            end
        end  
        
        
        
        
        function [GRAD] = POSTgrad2D( obj , u , NI , uP0x , uP0y )
            %POSTgrad2D: Calculates the gradient of a scalar u field in surface domain NI.
            %       GRAD = (uP0x*du/dx,uP0y*du/dy) 
            %INPUT
            %	- u: scalar field P0. 
            %   - NI: integer. Surface domain identifier. 
            %	- uP0x: a P0 scalar field. If it is a number, it is kept constant acros domain NI. Default 1. 
            %	- uP0y: a P0 scalar field. If it is a number, it is kept constant acros domain NI. Default 1	
            %OUTPUT: 
            %	- GRAD: a bidimensional vector field with the gradient. 

                tic;
                if ~exist('uP0x','var')
                    uP0x = 1;
                end   
                if ~exist('uP0y','var')
                    uP0y = 1;
                end   
                GRADx = uP0x.*( obj.create2DDerivativeP0( NI , 'dx' )*u );
                GRADy = uP0y.*( obj.create2DDerivativeP0( NI , 'dy' )*u );
                GRAD = [ GRADx ; GRADy ];    
                CPU = toc;    
                if obj.verbosity == 1
                    fprintf( 'Grad evaluated in %f s.\n' , CPU );
                end
        end
        
        
        
        
        function [CURL] = POSTcurl2D( obj , uz , NI , uP0x , uP0y )
            %POSTcurl2D: Calculates the curl of a vector u field in a surface domain NI within the plane XY. The curl operator is applied to a vector perpendicular to XY plane: (0,0,uz).
            %               CURL = ( uP0x*duz/dy ; -uP0y*duz/dx) );
            %INPUT
            %	- u: vector field z component (0,0,u). 
            %   - NI: integer. Surface domain identifier. 
            %	- uP0x: a P0 scalar field. If it is a number, it is kept constant acros domain NI. Default 1. 
            %	- uP0y: a P0 scalar field. If it is a number, it is kept constant acros domain NI. Default 1	
            %OUTPUT: 
            %	- CURL: a bidimensional vector field with the gradient. 
            
            tic;
            if ~exist('uP0x','var')
                uP0x = 1;
            end

            if ~exist('uP0y','var')
                uP0y = uP0x;
            end

            %para evaluar el rotacional en los nodos
            Mdx = obj.create2DDerivativeP0( NI , 'dx' );
            Mdy = obj.create2DDerivativeP0( NI , 'dy' );

            CURL = [ uP0x.*(Mdy*uz) ; uP0y.*(-Mdx*uz) ];

            CPU = toc;
            if obj.verbosity == 1
                fprintf( 'Curl 2D evaluated in %f s.\n' , CPU );
            end
        end

        
        
        
        function uP1 = POSTprojectP0toP1( obj , uP0 , NI , Ddim )
            %POSTprojectP0toP1(uP0 , NI , Ddim ): Interpolates (projects) the field P0, defined in the elements of domain NI, to a P1 field defined in the nodes of domain NI. 
                %INPUT
                %	- uP0: scalar or vector field defined in the domain elements (P0). 
                %   - NI: integer. Domain identifier. 
                %	- Ddim: dimension of ID domain (1 or ‘Line’ for curves, 2 or ‘Surface’ for surfaces). 
                %OUTPUT: 
                %	- uP1: interpolated (projected) scalar or vector field defined in the domain nodes (P1).

            switch Ddim
                case {2,'SurfaceIntegral'}
                    dimensions = length( uP0 )/obj.msh.nbTRI;
                case {1,'LineIntegral'}
                    dimensions = length( uP0 )/obj.msh.nbLINE;
            end
            
            uP0 = reshape( uP0 , [] , dimensions );    
            
            M   = obj.createMatrix( Ddim , NI , 'v*u' );
            
            for iDim = 1:dimensions
                L = obj.createVector( Ddim , NI , uP0(:,iDim) );
                uP1(:,iDim) = obj.SOLVER( M , L );
            end            

            uP1 = reshape( uP1 , [] , 1 );
            
        end 
        
        
        
        function [mean1,sum1,size1] = POSTmean( obj , DIM , NID , u )                
                     
  
            
               %u is P0 or P1
               if size( u , 1 )==obj.msh.nbNod
                    sum1 = obj.createVector( DIM , NID ).'*u;
                    size1 = sum( obj.createVector( DIM , NID ) );
                    mean1 = full( sum1 / size1 );

               else %if is P0
                    sum1 = sum( obj.createVector( DIM , NID , u ) );
                    size1 = sum( obj.createVector( DIM , NID ) );
                    mean1 = full( sum1 / size1 );
               end
        end
        

        
        
        
        
        

        

        %PLOT********************************************************

        function POSTplotGMSH2Dscalar( obj , u , surfaceID , name )
        %  POSTplotGMSH2Dscalar: creates a .pos file to be read by gmsh and plot the u scalar field.  The u field can be a scalar P0 or P1 field. The field name is XDYScalarP0.pos, where X is field name and Y is the surface domain identifier. 
        %INPUT
        %	- u: Scalar field to plot.
        %           - SurfaceID: integer vector. Surface domain identifier where the field is going to be ploted. If 0, computed in all domains.
        %	- name: A name for the scalar field. e.g. temperature or electric potential. 


                %read the scale
                scala = obj.scale;  
                
                u = full( u );
                               
                %The coordinates of nodes
                X = obj.msh.POS1(:,1)/scala;
                Y = obj.msh.POS1(:,2)/scala;
                Z = obj.msh.POS1(:,3)/scala;

                %create output file name and open it
                outfile = sprintf( '%sD%dScalarP0.pos' , name , surfaceID(1) );
                fileID = fopen(outfile,'w');  

                %Elements in domain NI
                elementos = MSHgetElementsNodes( obj , 2 , surfaceID );

                %load the triangles in the mesh
                TRI = obj.msh.TRIANGLES;
                
                %each node
                n1 = TRI(elementos,1);
                n2 = TRI(elementos,2);
                n3 = TRI(elementos,3);
                
                
                %coordinates of each triangle
                x1 = X(TRI(elementos,1));
                x2 = X(TRI(elementos,2));
                x3 = X(TRI(elementos,3));
                y1 = Y(TRI(elementos,1));
                y2 = Y(TRI(elementos,2));
                y3 = Y(TRI(elementos,3));
                z1 = Z(TRI(elementos,1));
                z2 = Z(TRI(elementos,2));
                z3 = Z(TRI(elementos,3));

                %u is P1, otherwhise (else) P0
                if length( u ) == obj.msh.nbNod
                        %P1
                        u1 = u( n1 );
                        u2 = u( n2 );
                        u3 = u( n3 );
                else%P0
                       u1 = u( elementos );
                       u2 = u1;
                       u3 = u1;
                end

                M = [x1 , y1 , z1 , x2 , y2 , z2 , x3 , y3 , z3  , u1 , u2 , u3 ]';

                %write
                fprintf(fileID,'View "%s 2D scalar field P0" {\n',name);
                fprintf(fileID,'ST(%e, %e, %e, %e, %e, %e, %e, %e, %e ){%e, %e, %e };\n', M );
                fprintf(fileID,'};');                
                
                %close
                fclose(fileID);

            end
            
            
            function POSTplotGMSH2Dvector( obj , u , surfaceID , name )
            %POSTplotGMSH2Dvector: creates a .pos file to be read by gmsh and plot the u vector field.  The u field can be a vector P0 or P1 field. The field name is XDYvectorP0.pos, where X is field name and Y is the surface domain identifier. 
            %INPUT
            %	- u: Vector field to plot.
            %           - SurfaceID: integer vector. Surface domain identifier where the field is going to be ploted. If 0, computed in all domains.
            %	- name: A name for the vector field. e.g. magnetic field.

                %read the scale
                scala = obj.scale;              

                %ad a null 3rd dimension
                u = [u ; zeros( length(u)/2 ,1 ) ]; 
                               
                %las coordenadas
                X = obj.msh.POS1(:,1)/scala;
                Y = obj.msh.POS1(:,2)/scala;
                Z = obj.msh.POS1(:,3)/scala;

                %create and open file .pos
                outfile = sprintf( '%sD%dvectorP0.pos' , name , surfaceID(1) );
                fileID = fopen(outfile,'w');
                
                                
                %the u components
                    u = reshape( u , [] , 3 );
                    ux = u( : , 1 );
                    uy = u( : , 2 );   
                    uz = u( : , 3 );  

                %Elements in domain NI
                elementos = MSHgetElementsNodes( obj , 2 , surfaceID(1) );

                %cargamos los tetaedros de msh
                TRI = obj.msh.TRIANGLES;
                
                %the index of each node
                n1 = TRI(elementos,1);
                n2 = TRI(elementos,2);
                n3 = TRI(elementos,3);
                
                %for each point in the triangle the coordinates
                x1 = X(n1);x2 = X(n2);x3 = X(n3);
                y1 = Y(n1);y2 = Y(n2);y3 = Y(n3);
                z1 = Z(n1);z2 = Z(n2);z3 = Z(n3);

               %u is P1, otherwhise (else) P0
               if size( u , 1 )==obj.msh.nbNod %P1
                   ux1 = ux( n1 );ux2 = ux( n2 );ux3 = ux( n3 );
                   uy1 = uy( n1 );uy2 = uy( n2 );uy3 = uy( n3 );
                   uz1 = uz( n1 );uz2 = uz( n2 );uz3 = uz( n3 );
               else %P0
                   ux1 = ux( elementos );ux2 = ux1; ux3 = ux1;
                   uy1 = uy( elementos );uy2 = uy1; uy3 = uy1;
                   uz1 = uz( elementos );uz2 = uz1; uz3 = uz1;
               end
               
               M = [x1 , y1 , z1 , x2 , y2 , z2 , x3 , y3 , z3  , ux1 , uy1 , uz1 , ux2 , uy2 , uz2 , ux3 , uy3 , uz3 ]';

                %write file
                fprintf(fileID,'View "%s 2D vector field P0" {\n',name);
                fprintf(fileID,'VT(%e, %e, %e, %e, %e, %e, %e, %e, %e ){%e, %e, %e, %e, %e, %e, %e, %e, %e };\n', M );

                %close file
                fprintf(fileID,'};');
                fclose(fileID);
            end
            
                      
            
            
            function POSTplotScalar2D( obj , u , surfaceID )
                %POSTplotScalar2D: plots the scalar field u into the surface domain
                %identified by surfaceID: 
                %INPUT: 
                %   - u: scalar field
                %   - surfaceID: integer. If 0, plots for all domains.
               
                 %points in surfaceID
                [ elementos , ~ ] = MSHgetElementsNodes( obj , 2 , surfaceID );
                
                %u is P0 or P1
               if size( u , 1 )==obj.msh.nbNod
                    u1 = u( obj.msh.TRIANGLES(elementos,1) );
                    u2 = u( obj.msh.TRIANGLES(elementos,2) );
                    u3 = u( obj.msh.TRIANGLES(elementos,3) );
                    U = [u1,u2,u3]';
               else %if is P0
                    u1 = u( elementos );
                    u2 = u( elementos );
                    u3 = u( elementos );
                    U = [u1,u2,u3]';
               end
                                
                x1 = obj.msh.POS1(obj.msh.TRIANGLES(elementos,1),1);
                x2 = obj.msh.POS1(obj.msh.TRIANGLES(elementos,2),1);
                x3 = obj.msh.POS1(obj.msh.TRIANGLES(elementos,3),1);
                y1 = obj.msh.POS1(obj.msh.TRIANGLES(elementos,1),2);
                y2 = obj.msh.POS1(obj.msh.TRIANGLES(elementos,2),2);
                y3 = obj.msh.POS1(obj.msh.TRIANGLES(elementos,3),2);
                X = [x1,x2,x3]';
                Y = [y1,y2,y3]';
                
                patch( X , Y , U , 'EdgeColor' , 'none');
                axis equal;        
            end
            
            
            
            
            
            
            
            function y = POSTplotScalar1D( obj , u , pointID )
                %POSTplotScalar1D: plots the scalar at point indentified by pointID: 
                %INPUT: 
                %   - u: scalar field
                %   - surfaceID: integer. If 0, plots for all domains.
               
                 %points in surfaceID
                [ ~ , node ] = MSHgetElementsNodes( obj , 0 , pointID );
                
                y = u( node );                     
            end





            function POSTplotVector2D( obj , u , surfaceID , geoFile )
                %POSTplotVector2D: plots the bidimensional vector field u into the surface domain
                %identified by surfaceID: 
                %INPUT: 
                %   - u: scalar field
                %   - surfaceID: integer. If 0, plots for all domains.
                %   - sampling: >1 integer. Plots the arrow each sampling element. Default 1, i.e. no sampling. 
                
                if ~exist( 'geoFile' , 'var' )
                    sampling = 1;    
                end
                
                vFEM = femSPACE2D( geoFile , 1 , obj.scale);                
                %points in surfaceID
                [ elementos , nodos ] = MSHgetElementsNodes( obj , 2 , surfaceID );                
                U = reshape( u , [] , 2 );
                %u is P0 or P1
               if size( U , 1 )==obj.msh.nbNod
                    x = obj.msh.POS1(nodos,1);
                    y = obj.msh.POS1(nodos,2);
                    ux = U( nodos , 1 );
                    uy = U( nodos , 2 );
                    xq = vFEM.msh.POS1(:,1);
                    yq = vFEM.msh.POS1(:,2);

                    
               else %if is P0
                    x = obj.msh.D2.rc(elementos,1);
                    y = obj.msh.D2.rc(elementos,2);
                    ux = U( elementos , 1 );
                    uy = U( elementos , 2 );
                    xq = vFEM.msh.D2.rc(:,1);
                    yq = vFEM.msh.D2.rc(:,2);
               end
               uxq = griddata(x,y,ux,xq,yq );
               uyq = griddata(x,y,uy,xq,yq );
               
               quiver( xq , yq , uxq , uyq , 'r' , 'LineWidth' , 2); 
               %axis equal;
            end
            
            
            
            
            
           function [v,s1] = POSTescalarLinea2D( obj , u , r1 , r2 )
           % POSTescalarLinea2D: Plotea el campo escalar u  en una linea
           % con 100 puntos. La linea va de P1 a P2.
               %Node coordinates
               x = obj.msh.POS1(:,1);
               y = obj.msh.POS1(:,2);

               %Teh MC of the triangles elements
               xc = obj.msh.D2.rc(:,1);
               yc = obj.msh.D2.rc(:,2);

               %points in the curve
               s = linspace(0,1,100);

               %distance
               dist = sqrt( sum( (r2-r1).^2 ) );

               %evaluate the field in
               xq = s*(r2(1)-r1(1))+r1(1);
               yq = s*(r2(2)-r1(2))+r1(2);

               %u is P1
               if length( u )==obj.msh.nbNod
                   v = griddata(x,y,u,xq,yq);
               else %u is PO
                   v = griddata(xc,yc,u,xq,yq);
               end

               s1 = s*dist;

           end
        
           %FUNCTIONS TO MANIPULATE THE u FIELDS (POST)********************************        
        
        function ui = FIELD_component( obj , uvector , component )
            %FIELD_component ( uvector , component): Returns u, a scalar field which is the component of uvector, a vectorial field. 
            %INPUT
                %	- uvector: vector field. P0 or P1. 
                %           - component: character string. ‘x’ or ‘y’ 
                %OUTPUT: 
                % 	- ui: scalar field. Component or modulus of uvector. 

                u31 = reshape( uvector , [] , 2 );

                switch component
                    case {1,'x'}
                        ui = u31( : , 1 );
                    case {2,'y'}
                        ui = u31( : ,2 );  
                    case {0,'norm'} %vector modulus
                        ui = sqrt( obj.Field_Dot( uvector , uvector ) );
                end
        end
                
                
                
        function u = Field_Dot( obj , uA , uB )
            %FIELD_Dot2D: Returns u, a scalar field which is the scalar multiplication of aA and uB vectorial fields. Fields are P0 or P1. 
            %INPUT
                %	- uA: P0 or P1 vectorial field (2D)
                %   - uB: P0 or P1 vectorial field (2D)
                %OUTPUT: 
                % 	- u: scalar P0 or P1 field. 


                uA22 = reshape( uA , [] , 2 );
                uB22 = reshape( uB , [] , 2 );

                uAx = uA22(:,1);
                uAy = uA22(:,2);

                uBx = uB22(:,1);
                uBy = uB22(:,2);

                %resultado
                u = uAx.*uBx + uAy.*uBy;
        end      
        
        
        
        
              
        
                
        
                
        
        
        function u = FIELD_Cross( uA , uB )
            %FIELD_Cross: Returns u, a scalar field which is the cross product of aA and uB vectorial fields. Fields are P0 or P1. 
            %INPUT
                %	- uA: P0 or P1 vectorial field (2D)
                %   - uB: P0 or P1 vectorial field (2D)
                %OUTPUT: 
                % 	- u: P0 or P1 vectorial field (2D)


                uA33 = reshape( uA , [] , 2 );
                uB33 = reshape( uB , [] , 2 );

                uAx = uA33(:,1);
                uAy = uA33(:,2);

                uBx = uB33(:,1);
                uBy = uB33(:,2);

                %resultado
                u = uAx.*uBy - uAy.*uBx; 
        end
        
        
        
        function u = FIELD_constantP0( obj , dim , NI , value )
            %FIELD_constantP0: Returns u, a constant vectorial or scalar P0 field defined in domain NI. 
            %INPUT
                %	- obj: femSPACE2D object.
                %   - dim: dimensions of NI domain. If it is a curve (1D), 1 or ‘Line’. If it is a surface, 2 or ‘Surface' 
                %   - NI: integer. Domain identifier. If 0, u is defined in all domains. 
                %   - value: integer or line vector: Constant value of u field
                %OUTPUT: 
                % 	- u: scalar or vectorial P0 field. 

                %NI domain elements indexs
                [ ielements , ~ ] = MSHgetElementsNodes( obj , dim , NI );

                %dimensions
                switch dim

                    case {1,'Line'}

                         u1 = zeros( obj.msh.nbLINE , 1 );

                    case {2,'Surface'}

                         u1 = zeros( obj.msh.nbTRI , 1 );            
                end

                %elements in NI are 1.
                u1( ielements ) = 1;

                %write the value.
                u = kron( value' , u1 );
        end 
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        function u = FIELD_constantP1( obj , dim , NI , valor )
        %FIELD_constantP1: Returns u, a constant vectorial or scalar P1 field defined in domain NI. 
        %INPUT
            %	- obj: femSPACE2D object.
            %           - dim: dimensions of NI domain. If it is a curve (1D), 1 or ‘Line’. If it is a surface, 2 or ‘Surface' 
        %- value: integer or line vector: Constant value of u field
        % NI: integer. Domain identifier. If 0, u is defined in all domains. 
            %OUTPUT: 
            % 	- u: scalar or vectorial P1 field. 

            %NI domain elements indexs
            [ ~ , inodos ] = MSHgetElementsNodes( obj , dim , NI );

            u1 = zeros( obj.msh.nbNod , 1 );

            %%elements in NI are 1.
            u1( inodos ) = 1;

            %write the values
            u = kron( valor' , u1 );
        end
            
    end
    
end

