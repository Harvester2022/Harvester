classdef MagnetostrictiveHarvester2D <handle
    %MagnetostrictiveHarvester2D Summary of this class goes here
    
    properties
        %fem space for elasticyty
        femSPACE_elasticity;
        %fem space for electromagnetic
        femSPACE_em;
        
        %to save the scale and gmsh file
        scale = 1;
        gmshFile = "";
        
        %Properties
        Prop_Elasticity; 
        Prop_Econdictivity;
        Prop_Permeability;
        Prop_Magnetostrictive;
        Prop_NoLinearPermeability;
        
        Geometry;
        
        %width of the system (OZ)
        Lz = 5e-3;
        
        %for the coils
        Coil = struct( 'WireSection' , 0.7854e-6 , 'Nturn' , 10 , 'sigma' , 60e6 );
        
        %For timing
        tvec;
        Periods = 2;             
        omega = 2*pi*60;
        Steps = 200;

        %for the magnetostrictive materials
        Magnets;
        
        %for the mass tips
        MassTips;
        
        %Air domain
        id_material_air = 1;
        
        %infinite boundary id (LINE)
        id_infty = 1;
        
        %Static line id:
        id_fix = 1;
        
        %physical constants
        mu0 = 4*pi*1e-7;
        g   = 9.81;
        RL = 1000e3;%Load resistance
        az = 1; %aceleration amplitude (m/s2)
        %damping rayleigh. Alpha (s^-1). Beta (s)
        alpha = 0;
        beta  = 0;
        
        %P0 fields
        P0FIELDS;
        
        %P1 fields
        P1FIELDS;
        
        %FEM matrixs
        Matrix;
        
        %To store the solution
        Solution;
        
    end
    
    methods
        function obj = MagnetostrictiveHarvester2D( gmshFile , scale , Parameters )
        % MagnetostrictiveHarvester2D - Class constructor for creating a 2D Magnetostrictive Harvester object.
        %   This function initializes the object by loading a GMSH mesh file, scaling the geometry,
        %   and extracting geometric properties of different domains.
        %
        % Inputs:
        %   gmshFile   - String representing the path to the GMSH file.
        %   scale      - Scaling factor for the geometry.
        %   Parameters - (Optional) Structure containing additional parameters. If not provided, defaults are used.
        %
        % Outputs:
        %   obj        - Initialized object with loaded geometry and calculated properties.

        % Check if the Parameters input exists, if not, set default value
        if ~exist('Parameters','var')
            % Default parameter value
            Parameters.a = 1;
        end

        % Store the scale and gmshFile inputs into the object
        obj.scale = scale;
        obj.gmshFile = gmshFile;

        % Call the method to create the GMSH model based on the provided Parameters
        obj.createGMSH( Parameters );

        % Obtain a list of unique domains from the mesh file (triangular elements)
        DomainList = unique( obj.femSPACE_em.msh.TRIANGLES(:,4) );

        % Load geometrical properties for each domain in the DomainList
        for iD = 1:1:length( DomainList )
            % Retrieve the list of elements belonging to the current domain
            ielements = MSHgetElementsNodes( obj.femSPACE_em , 2 ,  DomainList(iD) );

            % Calculate and store the area of the current domain
            obj.Geometry{DomainList(iD)}.Area = full( sum( obj.femSPACE_em.createVector( 'SurfaceIntegral' , DomainList(iD) ) ) );

            % Compute the mean coordinates (centroid) of the elements within the domain
            rmz = full( mean( obj.femSPACE_em.msh.D2.rc(ielements,:) , 1 ) );

            % Store the centroid (x, y coordinates only) of the current domain
            obj.Geometry{DomainList(iD)}.MZ = rmz(:,1:2); 
        end
        end
        
        
        function createGMSH( obj , Parameters )
            % createGMSH - Initializes finite element spaces for elasticity and electromagnetic problems.
            %
            % This method creates two finite element space objects using the provided GMSH file:
            %   1. Elasticity Space (2D elements)
            %   2. Electromagnetic Space (1D elements)
            %
            % Inputs:
            %   obj        - The object instance.
            %   Parameters - Structure containing additional parameters for the femSPACE2D objects.

            % Create a 2D finite element space for elasticity (second argument = 2)
            obj.femSPACE_elasticity = femSPACE2D( obj.gmshFile , 2 , obj.scale , Parameters );

            % Create a 2D finite element space for electromagnetic analysis (second argument = 1)
            obj.femSPACE_em = femSPACE2D( obj.gmshFile , 1 , obj.scale , Parameters );

            % Set verbosity levels to zero (no printing of messages)
            obj.femSPACE_elasticity.verbosity = 0;
            obj.femSPACE_em.verbosity = 0;                        
        end
        
        
        
        function LoadHarvester2DProblem( obj )      
            % LoadHarvester2DProblem
            % This function initializes and loads material properties, fields, and other simulation parameters
            % for a 2D harvester problem, including elasticity, permeability, electrical conductivity, 
            % mass tips, and coil properties. It also sets up initial and guessed solutions for the simulation.
            %
            % Inputs:
            %   obj - The object that stores all necessary fields and parameters for the simulation
            %         (including material properties, geometries, etc.)

            % 0 field initialization for P0 and P1 fields
            obj.P0FIELDS.u0 = obj.femSPACE_em.FIELD_constantP0( 'Surface' , 0 , 0 );
            obj.P1FIELDS.u0 = obj.femSPACE_em.FIELD_constantP1( 'Surface' , 0 , 0 );

            % Scalar P0 fields with default material properties
            uYoung   = obj.femSPACE_em.FIELD_constantP0( 'Surface' , 0 , 1 );   % Young's modulus
            uPoisson = obj.femSPACE_em.FIELD_constantP0( 'Surface' , 0 , 0.0 ); % Poisson's ratio
            umux     = obj.femSPACE_em.FIELD_constantP0( 'Surface' , 0 , obj.mu0 ); % Permeability in x-direction
            umuy     = obj.femSPACE_em.FIELD_constantP0( 'Surface' , 0 , obj.mu0 ); % Permeability in y-direction
            urho     = obj.P0FIELDS.u0; % Density (initialized to zero field)
            usigma   = obj.P0FIELDS.u0; % Electrical conductivity (initialized to zero field)
            uMd      = obj.femSPACE_em.FIELD_constantP0( 'Surface' , 0 , [0,0] ); % Magnetization (initialized to zero field)

            % Timing setup - generate time vector for the simulation
            obj.tvec = linspace( 0 , obj.Periods*2*pi/obj.omega , obj.Steps*obj.Periods );

            % Loop to update material properties based on elasticity properties
            for i1 = 1:length( obj.Prop_Elasticity )
                % Adding elasticity properties (Young's modulus, Poisson's ratio, and density)
                uYoung = uYoung + obj.femSPACE_elasticity.FIELD_constantP0( 'Surface' , obj.Prop_Elasticity{i1}.id , obj.Prop_Elasticity{i1}.YoungModulus );   
                uPoisson = uPoisson + obj.femSPACE_elasticity.FIELD_constantP0( 'Surface' , obj.Prop_Elasticity{i1}.id , obj.Prop_Elasticity{i1}.PoissonCoefficient );
                urho = urho + obj.femSPACE_elasticity.FIELD_constantP0( 'Surface' , obj.Prop_Elasticity{i1}.id , obj.Prop_Elasticity{i1}.Density );
            end

            % Loop to update material properties based on permeability
            for i1 = 1:length( obj.Prop_Permeability )
                umux = umux + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.Prop_Permeability{i1}.id , obj.Prop_Permeability{i1}.murx*obj.mu0 - obj.mu0 );
                umuy = umuy + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.Prop_Permeability{i1}.id , obj.Prop_Permeability{i1}.mury*obj.mu0 - obj.mu0 );
            end

            % Loop to update material properties based on electrical conductivity
            for i1 = 1:length( obj.Prop_Econdictivity )
                usigma = usigma + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.Prop_Econdictivity{i1}.id , obj.Prop_Econdictivity{i1}.value );
            end

            % Loop to update magnetization field based on magnet properties
            for i1 = 1:length( obj.Magnets )
                uMd = uMd + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.Magnets{i1}.id , obj.Magnets{i1}.value );
            end

            % Handling density in mass tips
            for i1 = 1:length( obj.MassTips )
                % Delete the density in mass tips
                ieMassTip = MSHgetElementsNodes( obj.femSPACE_em , 'Surface' ,  obj.MassTips{i1}.id );
                urho( ieMassTip ) = 0;

                % Get section of the mass tips and calculate mass density
                S = MSHgetDomainSices( obj.femSPACE_elasticity , 'Surface' , obj.MassTips{i1}.id );
                rhoMT = obj.MassTips{i1}.value / S / obj.Lz; % mass tip density
                urho = urho + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.MassTips{i1}.id , rhoMT );                
            end

            % Handling coils: Updating conductivity properties for coils
            obj.Coil = CREATEcoil2D( obj , obj.Coil );
            usigma = usigma + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.Coil.In.id   , obj.Coil.In.sigmaEFF );
            usigma = usigma + obj.femSPACE_em.FIELD_constantP0( 'Surface' , obj.Coil.Out.id  , obj.Coil.Out.sigmaEFF );

            % Magnetostrictive materials handling
            for i1 = 1:length( obj.Prop_Magnetostrictive )
                Hdata = [];%Reset for each material
                Sdata = [];%Reset for each material
                if isfield( obj.Prop_Magnetostrictive{i1} , 'id' )
                    % Read magnetostrictive data from file
                    Data = dlmread( obj.Prop_Magnetostrictive{i1}.datafile );
                    Hdata = Data(:,1); Hdata = [Hdata; 1000e3; 1000e3; 0; 0];
                    Sdata = -Data(:,2);Sdata = [Sdata; -100e6; 100e6; -100e6; 100e6];
                    Mdata = Data(:,3); Mdata = [Mdata; max(Mdata); max(Mdata); 0; 0];

                    % Interpolation of data for H and M
                    [HH, SS] = meshgrid( 0:5e3:1000e3 , -1e8:5e6:1e8 );
                    MM = griddata( Hdata , Sdata , Mdata , HH , SS );
                    MM(:,1) = 0; % Magnetization is null at H=0

                    % Calculate permeability based on magnetostrictive data
                    MUMU = obj.mu0*( 1 + MM./HH/obj.mu0 );
                    MUMU(:,1) = MUMU(:,2); % The susceptibility at H=0 (low field)

                    % Save calculated data
                    obj.Prop_Magnetostrictive{i1}.MM = MM;
                    obj.Prop_Magnetostrictive{i1}.HH = HH;
                    obj.Prop_Magnetostrictive{i1}.SS = SS;
                    obj.Prop_Magnetostrictive{i1}.MUMU = MUMU;

                    % Indices of triangle elements in the magnetostrictive material
                    obj.Prop_Magnetostrictive{i1}.ie = MSHgetElementsNodes( obj.femSPACE_em , 2 , obj.Prop_Magnetostrictive{i1}.id );
                end
            end

            % No-linear permeability handling            
            for i1 = 1:length( obj.Prop_NoLinearPermeability )
                BB = [];%Reset for each material
                if isfield( obj.Prop_NoLinearPermeability{i1} , 'id' )
                    % Read non-linear permeability data from file
                    Data = dlmread( obj.Prop_NoLinearPermeability{i1}.datafile );
                    Hdata = Data(:,1); % Magnetic field intensity
                    Bdata = Data(:,2); % Magnetization flux density (B)

                    % Interpolation of H and B data
                    HH = (0:1:10000)*1e3;
                    BB = min( interp1( Hdata , Bdata , HH , 'linear' , 'extrap' ) , Bdata(end) );
                    BB(1) = 0; % B is null at H=0
                    MUMU = BB ./ HH; % Permeability
                    MUMU(1) = MUMU(2); % The susceptibility at H=0 (low field)

                    % Save calculated data
                    obj.Prop_NoLinearPermeability{i1}.HH = HH;
                    obj.Prop_NoLinearPermeability{i1}.BB = BB;
                    obj.Prop_NoLinearPermeability{i1}.MUMU = MUMU;

                    % Indices of triangle elements in the non-linear permeability region
                    obj.Prop_NoLinearPermeability{i1}.ie = MSHgetElementsNodes( obj.femSPACE_em , 2 , obj.Prop_NoLinearPermeability{i1}.id );
                end
            end

            % Calculate Cijkl tensor (elastic strain tensor)
            uCijkl = ENSEMBLE_strain_C_tensor( 2 , uYoung , uPoisson );

            % Save all material properties and fields for later use
            obj.P0FIELDS.uYoung   = uYoung;
            obj.P0FIELDS.uPoisson = uPoisson;
            obj.P0FIELDS.urho     = urho;
            obj.P0FIELDS.umux     = umux;
            obj.P0FIELDS.umuy     = umuy;
            obj.P0FIELDS.uEconductivity = usigma;
            obj.P0FIELDS.uMd            = uMd;
            obj.P0FIELDS.uCijkl         = uCijkl;

            % Assemble matrices for later use (e.g., for solving the system)
            obj.ensebbleMatrix();

            % Remove initial and guessed solutions to reset the simulation
            obj.Solution = [];

            % Initialize static and dynamic solution fields
            obj.Solution.StaticMagentic.uH = [ obj.P0FIELDS.u0 ; obj.P0FIELDS.u0 ];
            obj.Solution.StaticMagentic.umux = obj.P0FIELDS.u0;
            obj.Solution.StaticMagentic.umuy = obj.P0FIELDS.u0;
            obj.Solution.Dynamic.uAT(:,1) = [ obj.P1FIELDS.u0 ;0;0;0];
            obj.Solution.Dynamic.uB(:,1) = [ obj.P0FIELDS.u0 ; obj.P0FIELDS.u0 ];
            obj.Solution.Dynamic.uH(:,1) = [ obj.P0FIELDS.u0 ; obj.P0FIELDS.u0 ];
            obj.Solution.Dynamic.ux(:,1) = obj.P1FIELDS.u0;
            obj.Solution.Dynamic.uy(:,1) = obj.P1FIELDS.u0;
            obj.Solution.Dynamic.vx(:,1) = obj.P1FIELDS.u0;
            obj.Solution.Dynamic.vy(:,1) = obj.P1FIELDS.u0;
            obj.Solution.Dynamic.umux(:,1) = obj.P0FIELDS.umux;
            obj.Solution.Dynamic.umuy(:,1) = obj.P0FIELDS.umuy;
        end   
        
        
        function ensebbleMatrix( obj )         
            % ensebbleMatrix
            % This function assembles matrices used for solving the system later on.
            % It includes matrices related to elasticity, mass, damping, magnetic 
            % loading, and non-linear elements. The assembled matrices are saved 
            % as properties of the object for later use in the simulation.
            %
            % Inputs:
            %   obj - The object that stores all the necessary fields, such as 
            %         material properties, geometries, and other parameters.

            % Elasticity - Assembling the stiffness matrix
            obj.Matrix.MS = +obj.femSPACE_elasticity.ENSEMBLE_Stiffnes2D( 0 , obj.P0FIELDS.uCijkl );

            % Elasticity - Assembling the load vector (for surface integral)
            obj.Matrix.Lg = +obj.femSPACE_elasticity.ENSEMBLE_Load( 'SurfaceIntegral' , 0 , kron( [0;-1] , obj.P0FIELDS.urho*obj.g ) );

            % Elasticity - Assembling the mass matrix (for surface integral)
            obj.Matrix.Elastic.Mmass    = obj.femSPACE_elasticity.ENSEMBLE_Mass2D( 'SurfaceIntegral' , 0 , obj.P0FIELDS.urho );

            % Elasticity - Assembling the damping matrix, which is a function of 
            % the mass matrix and a damping coefficient (alpha and beta parameters)
            obj.Matrix.Elastic.Mdamping = obj.Matrix.Elastic.Mmass*obj.alpha - obj.Matrix.MS*obj.beta;

            % Magnetic - Assembling the static magnetic load based on the magnetization field
            obj.Matrix.Magnetic.B_Jmag  = CREATE_Js_Mag2D( obj.femSPACE_em , 0 , obj.P0FIELDS.uMd );

            % Non-linear elements handling: Collecting indices for non-linear 
            % materials (magnetostrictive and non-linear permeability)
            obj.Matrix.ie_nl = [];

            % Collect indices for magnetostrictive materials
            for i1 = 1:length( obj.Prop_Magnetostrictive )                
                obj.Matrix.ie_nl = [obj.Matrix.ie_nl ; obj.Prop_Magnetostrictive{i1}.ie ];
            end

            % Collect indices for non-linear permeability materials
            for i1 = 1:length( obj.Prop_NoLinearPermeability )                
                obj.Matrix.ie_nl = [obj.Matrix.ie_nl ; obj.Prop_NoLinearPermeability{i1}.ie ] ;
            end

            % Flags - Check if there are any non-linear magnetic materials (e.g., 
            % magnetostrictive or non-linear permeability) and set the appropriate flag
            if length( obj.Prop_Magnetostrictive ) + length( obj.Prop_NoLinearPermeability ) > 0
                obj.Matrix.Flags.NoLinearMagnetic = 1;  % Non-linear magnetic materials exist
            else
                obj.Matrix.Flags.NoLinearMagnetic = 0;  % No non-linear magnetic materials
            end
        end

        
        
        function Solve_Estatic_Elasticy( obj )
            % Solve_Estatic_Elasticity
            % This function solves the static elasticity problem for the system.
            % It assembles the system of equations, applies the necessary boundary 
            % conditions, and computes the displacement and stress fields.
            %
            % Inputs:
            %   obj - The object that contains the necessary matrices, material properties, 
            %         boundary conditions, and fields for solving the static elasticity problem.

            % Assemble the system of equations
            A = -obj.Matrix.MS;   % Stiffness matrix (negative sign for consistency with equilibrium)
            B = +obj.Matrix.Lg;   % Load vector (positive sign as it's the load term)

            % Apply Dirichlet boundary conditions (null displacement on fixed line)
            % 'Line' specifies that the condition is applied along a line, and the 
            % displacements on the line [obj.id_fix, obj.id_infty] are set to zero.
            [A, B] = obj.femSPACE_elasticity.SOLVERapplyDirichlet( A, B, 'Line', [obj.id_fix, obj.id_infty], [0, 0] );

            % Solve the system of equations using the solver provided in the femSPACE_elasticity object
            u = obj.femSPACE_elasticity.SOLVER( A, B );

            % Extract the displacement components in x and y directions
            ux = obj.femSPACE_elasticity.FIELD_component( u, 'x' );
            uy = obj.femSPACE_elasticity.FIELD_component( u, 'y' );

            % Post-processing: Compute the stress (sigma) fields based on the displacement field
            uCijkl = obj.P0FIELDS.uCijkl;  % Material stiffness tensor
            uP0Sigma = obj.femSPACE_elasticity.POST_2DSIGMA_P0( u, 0, uCijkl );

            % Save the results in the solution structure
            obj.Solution.StaticElasticity.ux = ux;  
            obj.Solution.StaticElasticity.uy = uy; 
            obj.Solution.StaticElasticity.usigmaXX = obj.femSPACE_elasticity.FIELD_component( uP0Sigma(:, 1), 'x' );
            obj.Solution.StaticElasticity.usigmaXY = obj.femSPACE_elasticity.FIELD_component( uP0Sigma(:, 1), 'y' );
            obj.Solution.StaticElasticity.usigmaYY = obj.femSPACE_elasticity.FIELD_component( uP0Sigma(:, 2), 'y' );

            % Print a message
            fprintf( 'Static Elastic Fields Solved.\n' );
        end

        
        function Solve_Harmonic_Elasticyty( obj )
            % Solve_Harmonic_Elasticity
            % This function solves the harmonic elasticity problem for the system.
            % It computes the displacement, stress fields, and updates the mesh for 
            % the harmonic case considering damping, mass, and stiffness.
            %
            % This function requires that the `Solve_Estatic_Elasticy` function has 
            % been executed previously to solve for the static elasticity fields. 
            % The static elasticity results are used as a base for the harmonic solution.
            %
            % Inputs:
            %   obj - The object that contains the necessary matrices, material properties, 
            %         boundary conditions, and fields for solving the harmonic elasticity problem.
            %
            % Outputs:
            %   obj.Solution.Dynamic - Contains the dynamic displacement and stress fields, 
            %                           including harmonic components.
            %   obj.Solution.HarmonicElasticity - Stores the harmonic elasticity solution 
            %                                      (displacements and stresses).
            %   obj.femSPACE_em.msh.POS1 - The updated mesh after applying the displacements.

            % Assemble the load matrix
            MS = +obj.Matrix.MS;  % Stiffness matrix
            Lg = +obj.Matrix.Lg;  % Load vector
            % Assemble the dynamic system matrix (including mass, stiffness, and damping)
            Aha = -MS - obj.omega^2 * obj.Matrix.Elastic.Mmass + 1i * obj.omega * obj.Matrix.Elastic.Mdamping;
            Bha = +Lg * obj.az;   % Right-hand side (load) vector

            % Apply Dirichlet boundary conditions (null displacement on the fixed line)
            % 'Line' specifies that the condition is applied along a line, and the 
            % displacements on the line [obj.id_fix, obj.id_infty] are set to zero.
            [A, B] = obj.femSPACE_elasticity.SOLVERapplyDirichlet( Aha, Bha, 'Line', [obj.id_fix, obj.id_infty], [0, 0] );

            % Solve the system of equations for the displacements
            u = obj.femSPACE_elasticity.SOLVER( A, B );

            % Extract displacement components in x and y directions
            ux = obj.femSPACE_elasticity.FIELD_component( u, 'x' );
            uy = obj.femSPACE_elasticity.FIELD_component( u, 'y' );

            % Post-processing: Compute the stress components based on the displacement field
            uCijkl = obj.P0FIELDS.uCijkl;  % Material stiffness tensor
            uP0Sigma = obj.femSPACE_elasticity.POST_2DSIGMA_P0( u, 0, uCijkl );
            uSxx = obj.femSPACE_elasticity.FIELD_component( uP0Sigma(:, 1), 'x' );
            uSxy = obj.femSPACE_elasticity.FIELD_component( uP0Sigma(:, 1), 'y' );
            uSyy = obj.femSPACE_elasticity.FIELD_component( uP0Sigma(:, 2), 'y' );

            % Save results into the solution structure for dynamic analysis
            % Real part of displacement and stress fields for harmonic time dependence
            obj.Solution.Dynamic.ux = obj.Solution.StaticElasticity.ux + real(kron(exp(1i * obj.omega * obj.tvec), ux));   
            obj.Solution.Dynamic.uy = obj.Solution.StaticElasticity.uy + real(kron(exp(1i * obj.omega * obj.tvec), uy));
            obj.Solution.Dynamic.uSxx = obj.Solution.StaticElasticity.usigmaXX + real(kron(exp(1i * obj.omega * obj.tvec), uSxx));
            obj.Solution.Dynamic.uSxy = obj.Solution.StaticElasticity.usigmaXY + real(kron(exp(1i * obj.omega * obj.tvec), uSxy));
            obj.Solution.Dynamic.uSyy = obj.Solution.StaticElasticity.usigmaYY + real(kron(exp(1i * obj.omega * obj.tvec), uSyy));

            % Save harmonic elasticity results
            obj.Solution.HarmonicElasticity.u = u;
            obj.Solution.HarmonicElasticity.ux = ux;
            obj.Solution.HarmonicElasticity.uy = uy;
            obj.Solution.HarmonicElasticity.uSxx = uSxx;
            obj.Solution.HarmonicElasticity.uSxy = uSxy;
            obj.Solution.HarmonicElasticity.uSyy = uSyy;

            % Update the mesh to reflect the new displacement from the solution
            % Move the mesh nodes based on the computed displacements
            obj.femSPACE_em.msh.POS1(:, 1) = obj.femSPACE_em.msh.POS(:, 1) + obj.Solution.Dynamic.ux(:, 1);
            obj.femSPACE_em.msh.POS1(:, 2) = obj.femSPACE_em.msh.POS(:, 2) + obj.Solution.Dynamic.uy(:, 1);
            obj.femSPACE_em.prepareMesh();  % Prepare the mesh for next steps

            % Update matrices for later use (e.g., after mesh update)
            obj.ensebbleMatrix();

            % Print message indicating the harmonic elasticity problem has been solved
            fprintf('Harmonic Elastic Fields Solved.\n');
        end

        
        
        
        function Solve_VibrationModes( obj )
            % Solve_VibrationModes
            % This function solves for the vibration modes of the system, including 
            % the associated natural frequencies (in Hz). It computes the eigenmodes 
            % and vibration frequencies using the mass and stiffness matrices.
            %
            % Inputs:
            %   obj - The object containing the stiffness matrix (MS), mass matrix 
            %         (Mmass), boundary conditions (id_fix), and other relevant fields.
            %
            % Outputs:
            %   obj.Solution.VibrationModes.u - The vibration mode shapes (eigenvectors).
            %   obj.Solution.VibrationModes.vfreqs - The natural vibration frequencies 
            %                                         (in Hz) corresponding to the modes.

            % Load the stiffness matrix
            MS = +obj.Matrix.MS;  

            % Load the mass matrix
            Mss = obj.Matrix.Elastic.Mmass;  

            % Solve for the vibration modes (eigenvectors) and natural frequencies
            [v, u] = obj.femSPACE_elasticity.SOLVEmodes(MS, Mss, obj.id_fix); 

            % Compute the natural frequencies (in Hz) from the eigenvalues
            vfreq = 0.5 * (1i * obj.alpha - 1i * sqrt(4 * v + obj.alpha^2)) / (2 * pi);

            % Save results into the solution structure
            obj.Solution.VibrationModes.u = u;        % Vibration mode shapes (eigenvectors)
            obj.Solution.VibrationModes.vfreqs = vfreq;  % Natural vibration frequencies (Hz)

            % Print message indicating that vibration modes have been solved
            fprintf('Vibration modes Solved.\n');
        end

        
        
        
        function Solve_StaticMagnetic_Linear( obj )
            % Solve_StaticMagnetic_Linear
            % This function solves for the static magnetic field in a linear
            % system. It computes the magnetic vector potential (uA), magnetic flux density (uB),
            % and magnetic field (uH) based on the given conditions and boundary conditions. 
            % The solution is computed using finite element methods and post-processing to 
            % calculate relevant magnetic field quantities.
            %
            % Inputs:
            %   obj - The object containing matrices, boundary conditions, and necessary 
            %         fields related to magnetic properties for the computation.           %
            %
            % Notes:
            %   - If nonlinear magnetic and magnetostrictive materials are present, the nonlinear permeability (umux and umuy)
            %     is computed using the functions `noLinear_mux` and `noLinear_muy`.
            %   - The magnetic vector potential (uA) is calculated using a curl-curl formulation.
            %   - Dirichlet boundary conditions are applied at infinity for the null potential vector.
            %   - This function uses the finite element method (FEM) for solving the system.
            %   - The solution is calculated in both static and dynamic contexts.
            
            % Check for nonlinear magnetic materials
            if obj.Matrix.Flags.NoLinearMagnetic==1  
                % If nonlinear magnetic properties are present, compute the nonlinear permeability
                uH = obj.Solution.StaticMagentic.uH;   % Previously calculated magnetic field         
                uSxx = obj.Solution.Dynamic.uSxx(:,1); % Stress in the x-direction
                uSyy = obj.Solution.Dynamic.uSyy(:,1); % Stress in the y-direction           
                % Compute the nonlinear permeability based on the current magnetic field and stress components
                umux = obj.noLinear_mux( uH , uSxx );
                umuy = obj.noLinear_muy( uH , uSyy );
            else
                % If no nonlinear magnetic materials, use linear permeability
                umux = obj.Solution.Dynamic.umux(:,1);
                umuy = obj.Solution.Dynamic.umuy(:,1);  
            end
            
            % Construct the stiffness (CurlCurl) matrix for the 2D problem
            MRR = obj.femSPACE_em.ENSEMBLE_curlcurl2D( 0 , 1./umux , 1./umuy );
                       
            % Apply Dirichlet boundary conditions at infinity (null potential vector)
            [A,B] = obj.femSPACE_em.SOLVERapplyDirichlet( MRR , obj.Matrix.Magnetic.B_Jmag , 'Line' , obj.id_infty , 0 );
            
            % Solve for the magnetic vector potential (uA) z component
            uA = obj.femSPACE_em.SOLVER( A , B );
            
            %Save into solution
            obj.Solution.StaticMagentic.uA = uA; 
            % Post-process to calculate the magnetic flux density (uB) and magnetic field (uH)
            obj.Solution.StaticMagentic.uB = obj.femSPACE_em.POSTcurl2D( uA , 0 );
            obj.Solution.StaticMagentic.uH = obj.femSPACE_em.POSTcurl2D( uA , 0 , 1./umux , 1./umuy );
             % Save the computed nonlinear permeability values
            obj.Solution.StaticMagentic.umux = umux;
            obj.Solution.StaticMagentic.umuy = umuy;
            
            % Save dynamic solutions (used in later time-stepping or dynamic analysis)
            obj.Solution.Dynamic.uAT(:,1) = [ uA ; 0 ; 0 ; 0 ];
            obj.Solution.Dynamic.uB(:,1) = obj.femSPACE_em.POSTcurl2D( uA , 0 );
            obj.Solution.Dynamic.uH(:,1) = obj.femSPACE_em.POSTcurl2D( uA , 0 , 1./umux , 1./umuy );
            obj.Solution.Dynamic.umux(:,1) = umux;
            obj.Solution.Dynamic.umuy(:,1) = umuy;
            
            % Print message indicating the static magnetic field has been solved
            fprintf( 'Static magnetic fiels solved.\n' )
            
        end   
        
        
        
        
        function Solve_StaticMagnetic_noLinear( obj )
            % Solve_StaticMagnetic_noLinear: Iteratively solves for the static magnetic field
            % in a nonlinear system.
            %
            % This function updates the magnetic permeability (umux, umuy) iteratively by
            % solving the linearized magnetic field using `Solve_StaticMagnetic_Linear`. 
            % It continues updating until the relative error between the current and previous
            % magnetic permeability values is below a specified tolerance.
            %
            % The function works in the following steps:
            %   1. Initializes the tolerance (`tol`) and maximum iteration count (100).
            %   2. At each iteration, it computes the magnetic permeability values for the nonlinear
            %      elements based on the solution from the previous iteration.
            %   3. Calls `Solve_StaticMagnetic_Linear` to solve the linearized system.
            %   4. Computes the error between the current and previous permeability values.
            %   5. Stops iterating if the error is smaller than the tolerance.
            %   
            % Inputs:
            %   obj - The object that contains the parameters and fields for the magnetic solution.
            %
            % Outputs:
            %   The function updates the `obj.Solution.StaticMagentic` field with the computed
            %   magnetic field solutions and permeabilities.
            %
            % Example:
            %   obj.Solve_StaticMagnetic_noLinear();

            % Define tolerance for convergence
            tol = 1e-2;

            % Maximum number of iterations
            for i1 = 1:100
                % Store the current magnetic permeability values (before update)
                umu_k0 = [ obj.Solution.StaticMagentic.umux(  obj.Matrix.ie_nl ) ; obj.Solution.StaticMagentic.umuy( obj.Matrix.ie_nl ) ];

                % Solve for the static magnetic field with the current permeabilities (linearized)
                Solve_StaticMagnetic_Linear( obj );

                % Get the updated permeability values
                umu_k = [ obj.Solution.StaticMagentic.umux(  obj.Matrix.ie_nl ) ; obj.Solution.StaticMagentic.umuy( obj.Matrix.ie_nl ) ];

                % Compute the error between the previous and current permeability values
                error = sqrt( (umu_k0-umu_k).'*(umu_k0-umu_k) / (umu_k0.'*umu_k0) );

                % Display the iteration number and the current error
                fprintf( 'NonLinear it = %d. Error = %.9f\n' , i1 , error )

                % If the error is below the tolerance, stop iterating
                if error < tol
                    break
                end
            end
        end

        
        
        
        
        function Matgentic_time_step( obj ,  k )          
            % Matgentic_time_step: Updates the magnetic field and forces at each time step 
            % in a dynamic simulation, considering the magnetic material's nonlinear behavior.
            %
            % This function updates the magnetic field (`uAT`, `uB`, `uH`) at each time step 
            % by solving the system of equations that govern the evolution of the magnetic field.
            % The function also updates other necessary fields, such as the magnetic permeability 
            % (`umux`, `umuy`) and the magnetic force (`uFmag`), as well as the matrices needed for 
            % subsequent time steps. This method handles the nonlinear nature of the material 
            % through an iterative approach and by updating the material properties at each step.
            %
            % The main steps include:
            %   1. Retrieve the previous values of magnetic field components and stresses.
            %   2. Compute the magnetic permeability based on previous values.
            %   3. Compute the stiffness matrix for the current time step.
            %   4. Update coil information and the associated magnetic load.
            %   5. Assemble the system matrix `A` for the current time step.
            %   6. Solve the system of equations for the magnetic field and force.
            %   7. Store the computed values for use in future time steps.
            %   
            % Inputs:
            %   obj - The object that contains the fields and methods for solving the magnetic field.
            %   k   - The current time step index.
            %
            % Outputs:
            %   The function updates various fields within `obj.Solution.Dynamic`, including:
            %     - `uAT`: The vector potential solution.
            %     - `uB`: The magnetic flux density solution.
            %     - `uH`: The magnetic field strength solution.
            %     - `umux`, `umuy`: The magnetic permeability values.
            %     - `uM`: The magnetic field modification due to external influences.
            %     - `uFmag`: The magnetic force.

            % Previous step fields
            uSxx_0 = obj.Solution.Dynamic.uSxx(:,k-1);
            uSyy_0 = obj.Solution.Dynamic.uSyy(:,k-1);
            uH_0 = obj.Solution.Dynamic.uH(:,k-1);

            % Load mu. Evaluated with S(k-1) and H(k-1)
            umux_0 = obj.noLinear_mux( uH_0 , uSxx_0 );
            umuy_0 = obj.noLinear_muy( uH_0 , uSyy_0 );

            % If it is the first call (2nd step), obj.Matrix.MRR_0 does not exist
            if k==2
                % Move mesh (magnetic)
                umux_0 = obj.Solution.Dynamic.umux(:,1);
                umuy_0 = obj.Solution.Dynamic.umuy(:,1);
                obj.Matrix.MRR_0 = obj.femSPACE_em.ENSEMBLE_curlcurl2D( 0 , 1./umux_0 , 1./umuy_0 );
            end

            % Stiffness (CurlCurl) Matrix
            MRR   = obj.femSPACE_em.ENSEMBLE_curlcurl2D( 0 , 1./umux_0 , 1./umuy_0 );
            MRR_0 = obj.Matrix.MRR_0; % Previous stiffness matrix
            obj.Matrix.MRR_0 = MRR;

            % Mass matrix
            Mmass = obj.femSPACE_em.ENSEMBLE_Mass2D( 'SurfaceIntegral' , 0 , obj.P0FIELDS.uEconductivity );

            % Update coils
            Coil2 = CREATEcoil2D( obj , obj.Coil );

            % 1V currents
            BJ1Vin  = Coil2.In.BJs;
            BJ1Vout = Coil2.Out.BJs;
            EMFin   = Coil2.In.emfIntegral;
            EMFout  = Coil2.Out.emfIntegral;
            RdcIn   = Coil2.In.Rdc;
            RdcOut  = Coil2.Out.Rdc;

            % Magnetic Load
            BJmag   = CREATE_Js_Mag2D( obj.femSPACE_em , 0 , obj.P0FIELDS.uMd );
            u0P1 = zeros( size( BJ1Vin ) );

            % Time increment (regular)
            dt = obj.tvec(2) - obj.tvec(1);

            % Ensemble Matrix
            A  = [ +MRR  /2 + Mmass/dt, -BJ1Vin/2, -BJ1Vout/2, u0P1;
                   -EMFin./dt, -1/2, 0, -RdcIn/2;
                   -EMFout./dt, 0, -1/2, -RdcOut/2;
                   u0P1', +1/2, +1/2, -obj.RL/2 ];
            A0 = [ -MRR_0/2 + Mmass/dt, +BJ1Vin/2, +BJ1Vout/2, u0P1;
                   -EMFin./dt, +1/2, 0, +RdcIn/2;
                   -EMFout./dt, 0, +1/2, +RdcOut/2;
                   u0P1', -1/2, -1/2, +obj.RL/2 ];
            B0 = [ BJmag ; 0 ; 0 ; 0 ];

            % Compute the right-hand side vector
            B = A0 * obj.Solution.Dynamic.uAT(:,k-1) + B0;

            % Apply Dirichlet condition in the infinity (null potential vector)
            [A ,B] = obj.femSPACE_em.SOLVERapplyDirichlet( A, B, 'Line', obj.id_infty, 0 );      

            % Solve for the unknowns
            uAT = obj.femSPACE_em.SOLVER( A , B );
            uB = obj.femSPACE_em.POSTcurl2D( uAT(1:end-3) , 0 );
            uH = obj.femSPACE_em.POSTcurl2D( uAT(1:end-3) , 0 , 1./umux_0 , 1./umuy_0 );
            uM = uB/obj.mu0 - uH + obj.P0FIELDS.uMd;  

            % Store results
            obj.Solution.Dynamic.uAT(:,k) = uAT;
            obj.Solution.Dynamic.uB(:,k) = uB;
            obj.Solution.Dynamic.uH(:,k) = uH;
            obj.Solution.Dynamic.umux(:,k) = umux_0;
            obj.Solution.Dynamic.umuy(:,k) = umuy_0;                       
            obj.Solution.Dynamic.uM(:,k) = uM; 
        end       
        
        
        
        
        
        
        
         function Solve_Dynamic( obj , BooleanMoveMesh )
            % Solve_Dynamic: Solves the dynamic system over a range of time steps.
            % This function performs dynamic simulations by solving for the magnetic and
            % elastic forces at each time step, optionally moving the mesh based on the
            % displacement and updating the geometry.
            %
            % The function loops over the time vector, solving for the magnetic forces and
            % updating the mesh if needed. It calls the `Matgentic_time_step` method to
            % compute the magnetic force and other dynamic properties for each time step.
            %
            % The time steps are solved for all time points from 2 to the length of the time 
            % vector `tvec`. Optionally, the mesh can be moved based on the displacement 
            % calculated in previous steps.
            %
            % Inputs:
            %   obj               - The object containing the matrices, solution data, 
            %                        and methods for dynamic simulation.
            %   BooleanMoveMesh   - Optional flag to move the mesh. If set to 1 or 
            %                        'Move geometry', the mesh will be moved according
            %                        to the displacement. Default is 0 (no mesh movement).
            %
            % Outputs:
            %   The function does not return any values but updates the system state 
            %   based on the dynamic solution (i.e., updates the mesh and solution data).

            % Set default values if inputs are not provided
            if ~exist( 'BooleanMoveMesh' , 'var' )
                BooleanMoveMesh = 0;  % Default: do not move mesh
            end
            if ~exist( 'BooleanNoLinearElasticity' , 'var' )
                BooleanNoLinearElasticity = 0;  % Default: linear elasticity
            end

            % Loop over time steps
            for k = 2:length( obj.tvec )

                % Move mesh if the BooleanMoveMesh flag is set (either as 1 or 'Move geometry')
                if or( BooleanMoveMesh == 1 , BooleanMoveMesh == 'Move geometry' )
                    % Update mesh position with displacement values
                    obj.femSPACE_em.msh.POS1(:,1) = obj.femSPACE_em.msh.POS(:,1) + obj.Solution.Dynamic.ux(:,k-1);
                    obj.femSPACE_em.msh.POS1(:,2) = obj.femSPACE_em.msh.POS(:,2) + obj.Solution.Dynamic.uy(:,k-1);
                    obj.femSPACE_em.prepareMesh();   % Recalculate the mesh
                end

                % Solve for magnetic forces and update dynamics
                obj.Matgentic_time_step( k );

                % Print progress every 10% of the total time steps
                if mod(10*k, length( obj.tvec )) == 0
                    fprintf( 'Magnetic time step %d from %d solved.\n', k, length( obj.tvec ) );
                end            
            end

            % Final message indicating the completion of the dynamic solution
            fprintf( 'Dynamic time steps Solved\n' );
         end  
        
        

        
        
         function umux = noLinear_mux( obj , uH , uSxx )
            % noLinear_mux: Computes the nonlinear permeability in the x-direction.
            %
            % This function computes the nonlinear permeability `umux` in the x-direction based on 
            % the applied magnetic field `uH` and the strain `uSxx`. It considers the magnetostrictive 
            % properties and nonlinear permeability data, clamping the strain and field to specific maximum 
            % values before performing interpolation to get the permeability.
            %
            % Inputs:
            %   obj    - The object containing the magnetostrictive properties and permeability data.
            %   uH     - The magnetic field vector (2D).
            %   uSxx   - The strain in the x-direction.
            %
            % Outputs:
            %   umux   - The resulting nonlinear permeability in the x-direction.
            %
            % Process:
            %   - The strain `uSxx` and magnetic field `uHx` are clamped to defined maximum values.
            %   - Interpolation is performed on the magnetostrictive properties `HH`, `SS`, and `MUMU` to 
            %     calculate `umux` for each element based on the current strain and magnetic field.
            %   - The function also supports properties defined for nonlinear permeability, where 
            %     interpolation is used to calculate the permeability based on the magnetic field.

            % Initialize umux to the default permeability
            umux = obj.P0FIELDS.umux;            

            % Extract the x-component of the magnetic field
            uHx = obj.femSPACE_em.FIELD_component( uH , 'x' );           

            % Loop through the magnetostrictive properties and calculate the permeability
            for i1 = 1:length( obj.Prop_Magnetostrictive )                       
                HH = obj.Prop_Magnetostrictive{i1}.HH;
                SS = obj.Prop_Magnetostrictive{i1}.SS;
                MUMU = obj.Prop_Magnetostrictive{i1}.MUMU;
                ie = obj.Prop_Magnetostrictive{i1}.ie;

                % Clamp strain and magnetic field to predefined maximum values
                S_max = 1e8;
                H_max  = 1e5;
                uSxx(ie) = min(max( uSxx(ie) , -S_max ), S_max );
                uHx(ie)  = min(max( uHx(ie) , -H_max ), H_max );             

                % Interpolate to calculate umux
                umux(ie) = griddata( HH , SS , MUMU , abs( uHx(ie) ), uSxx(ie) );         
            end

            % For nonlinear permeability properties, perform interpolation
            for i1 = 1:length( obj.Prop_NoLinearPermeability )
                HH = obj.Prop_NoLinearPermeability{i1}.HH;
                MUMU = obj.Prop_NoLinearPermeability{i1}.MUMU;
                ie = obj.Prop_NoLinearPermeability{i1}.ie;
                umux(ie) = interp1( HH , MUMU , abs( uHx(ie) ) , 'linear' , MUMU(end) ); 
                %when HH is out of range, last value of MUMU is provided: MUMU(end)
            end
         end

        
                
         function umuy = noLinear_muy( obj , uH , uSyy )
            % noLinear_muy: same as noLinear_mux but in y

            % Initialize umuy to the default permeability
            umuy = obj.P0FIELDS.umuy;            

            % Extract the y-component of the magnetic field
            uHy = obj.femSPACE_em.FIELD_component( uH , 'y' );           

            % Loop through the magnetostrictive properties and calculate the permeability
            for i1 = 1:length( obj.Prop_Magnetostrictive )                       
                HH = obj.Prop_Magnetostrictive{i1}.HH;
                SS = obj.Prop_Magnetostrictive{i1}.SS;
                MUMU = obj.Prop_Magnetostrictive{i1}.MUMU;
                ie = obj.Prop_Magnetostrictive{i1}.ie;

                % Clamp strain and magnetic field to predefined maximum values
                S_max = 1e8;
                H_max  = 1e5;
                uSyy(ie) = min(max( uSyy(ie) , -S_max ), S_max );
                uHy(ie)  = min(max( uHy(ie) , -H_max ), H_max );             

                % Interpolate to calculate umuy
                umuy(ie) = griddata( HH , SS , MUMU , abs( uHy(ie) ), uSyy(ie) );            
            end

            % For nonlinear permeability properties, perform interpolation
            for i1 = 1:length( obj.Prop_NoLinearPermeability )
                HH = obj.Prop_NoLinearPermeability{i1}.HH;
                MUMU = obj.Prop_NoLinearPermeability{i1}.MUMU;
                ie = obj.Prop_NoLinearPermeability{i1}.ie;
                umuy(ie) = interp1( HH , MUMU , abs( uHy(ie) ) );             
            end
        end

        
        

        
        
        function Plot_Estatic_Elasticyty( obj , field , id_vec )
            % Plot_Estatic_Elasticyty: This function generates plots for static elasticity results.
            %
            % This function allows the visualization of different static elasticity quantities, such as
            % displacements and stress components, depending on the specified field. The output is a 2D plot
            % where the x-axis and y-axis correspond to the spatial coordinates, and the field of interest 
            % (e.g., displacement, stress) is visualized using a color scale.
            %
            % Inputs:
            %   obj      - The object that contains the solution data for the static elasticity problem.
            %   field    - A string specifying which field to plot. Possible values are:
            %              'u'    - Magnitude of displacement.
            %              'ux'   - Displacement in the x-direction.
            %              'uy'   - Displacement in the y-direction.
            %              'Sxx'  - Stress in the x-direction.
            %              'Sxy'  - Shear stress (xy-plane).
            %              'Syy'  - Stress in the y-direction.
            %   id_vec   - A vector of element IDs specifying the region or mesh area to be plotted.
            %
            % Outputs:
            %   - A 2D plot showing the selected field as a function of spatial coordinates (X and Y).
            %     The title of the plot will reflect the quantity being plotted, and the axes will be labeled 
            %     with units (if applicable).
            %
            % Process:
            %   - Based on the input `field`, the appropriate static elasticity result is selected from 
            %     the solution object (`obj.Solution.StaticElasticity`).
            %   - The result is then plotted using the `POSTplotScalar2D` function, which generates a 
            %     color map for the selected field.
            %
            % Example Usage:
            %   obj.Plot_Estatic_Elasticyty('u', id_vec);   % Plot displacement magnitude
            %   obj.Plot_Estatic_Elasticyty('Syy', id_vec); % Plot stress in y-direction

            % Switch based on the field input
            switch field
                case 'u'            
                    % Modulus of the displacement (magnitude)
                    u = obj.femSPACE_elasticity.FIELD_component( obj.Solution.StaticElasticity.u , 'norm' );
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Displacement $u$ (m)', 'Interpreter', 'latex')

                case 'ux'            
                    % X-displacement
                    obj.femSPACE_elasticity.POSTplotScalar2D( obj.Solution.StaticElasticity.ux , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Displacement $u_x$ (m)', 'Interpreter', 'latex')

                case 'uy'            
                    % Y-displacement
                    obj.femSPACE_elasticity.POSTplotScalar2D( obj.Solution.StaticElasticity.uy , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Displacement $u_y$ (m)', 'Interpreter', 'latex')

                case 'Sxx'    
                    % Stress in the x-direction (sigma_xx)
                    u = obj.Solution.StaticElasticity.usigmaXX;
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{xx}$ (Pa)', 'Interpreter', 'latex')                    

                case 'Sxy'    
                    % Shear stress (sigma_xy)
                    u = obj.Solution.StaticElasticity.usigmaXY;
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{xy}$ (Pa)', 'Interpreter', 'latex')

                case 'Syy'    
                    % Stress in the y-direction (sigma_yy)
                    u = obj.Solution.StaticElasticity.usigmaYY;
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{yy}$ (Pa)', 'Interpreter', 'latex')
            end
        end

            
            
        function Plot_Harmonic_Elasticyty( obj , field , id_vec )
            % Plot_Harmonic_Elasticyty: This function generates plots for harmonic elasticity results.
            %
            % This function allows the visualization of different harmonic elasticity quantities, such as
            % displacements and stress components, depending on the specified field. The output is a 2D plot
            % where the x-axis and y-axis correspond to the spatial coordinates, and the field of interest 
            % (e.g., displacement, stress) is visualized using a color scale.
            %
            % Inputs:
            %   obj      - The object that contains the solution data for the harmonic elasticity problem.
            %   field    - A string specifying which field to plot. Possible values are:
            %              'u'    - Magnitude of displacement.
            %              'ux'   - Displacement in the x-direction.
            %              'uy'   - Displacement in the y-direction.
            %              'Sxx'  - Stress in the x-direction.
            %              'Sxy'  - Shear stress (xy-plane).
            %              'Syy'  - Stress in the y-direction.
            %   id_vec   - A vector of element IDs specifying the region or mesh area to be plotted.
            %
            % Outputs:
            %   - A 2D plot showing the selected field as a function of spatial coordinates (X and Y).
            %     The title of the plot will reflect the quantity being plotted, and the axes will be labeled 
            %     with units (if applicable).
            %
            % Process:
            %   - Based on the input `field`, the appropriate harmonic elasticity result is selected from 
            %     the solution object (`obj.Solution.HarmonicElasticity`).


            % Switch based on the field input
            switch field
                case 'u'            
                    % Modulus of the displacement (magnitude)
                    u = obj.femSPACE_elasticity.FIELD_component( obj.Solution.HarmonicElasticity.u , 'norm' );
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Harmonic Displacement $u$ (m)', 'Interpreter', 'latex')

                case 'ux'            
                    % X-displacement
                    obj.femSPACE_elasticity.POSTplotScalar2D( obj.Solution.HarmonicElasticity.ux , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Harmonic Displacement $u_x$ (m)', 'Interpreter', 'latex')

                case 'uy'            
                    % Y-displacement
                    obj.femSPACE_elasticity.POSTplotScalar2D( abs( obj.Solution.HarmonicElasticity.uy ) , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Harmonic Displacement $u_y$ (m)', 'Interpreter', 'latex')

                case 'Sxx'    
                    % Stress in the x-direction (sigma_xx)
                    u = obj.Solution.HarmonicElasticity.uSxx;
                    obj.femSPACE_elasticity.POSTplotScalar2D( abs( u ) , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{xx}$ (Pa)', 'Interpreter', 'latex')                    

                case 'Sxy'    
                    % Shear stress (sigma_xy)
                    u = obj.Solution.HarmonicElasticity.uSxy;
                    obj.femSPACE_elasticity.POSTplotScalar2D( abs( u ), id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{xy}$ (Pa)', 'Interpreter', 'latex')

                case 'Syy'    
                    % Stress in the y-direction (sigma_yy)
                    u = obj.Solution.HarmonicElasticity.uSyy;
                    obj.femSPACE_elasticity.POSTplotScalar2D( abs( u ) , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{yy}$ (Pa)', 'Interpreter', 'latex')
            end
        end

        
              
         function y = Get_Harmonic_Elasticyty_Point( obj , field , id_vec )
            % Get_Harmonic_Elasticyty_Point: This function retrieves and plots harmonic elasticity quantities.
            %
            % This function provides a way to extract and plot different harmonic elasticity results for 
            % selected fields (such as displacement or stress components). Depending on the specified 
            % field, the function retrieves the corresponding data and generates a 1D or 2D plot for the 
            % given region of interest (defined by `id_vec`).
            %
            % Inputs:
            %   obj      - The object that contains the solution data for the harmonic elasticity problem.
            %   field    - A string specifying which field to retrieve and plot. Possible values are:
            %              'u'    - Magnitude of displacement.
            %              'ux'   - Displacement in the x-direction.
            %              'uy'   - Displacement in the y-direction.
            %              'Sxx'  - Stress in the x-direction.
            %              'Sxy'  - Shear stress (xy-plane).
            %              'Syy'  - Stress in the y-direction.
            %   id_vec   - A vector of element IDs specifying the region or mesh area to be plotted.
            %
            % Outputs:
            %   y        - The result of the field at the specified element IDs. This will vary based on the
            %              type of field. It could be a scalar or a matrix, depending on the field.
            %
            % Process:
            %   - Based on the `field` input, the function retrieves the corresponding harmonic elasticity 
            %     result from `obj.Solution.HarmonicElasticity` and plots the field.
            %
            % Example Usage:
            %   y = obj.Get_Harmonic_Elasticyty_Point('u', id_vec);   % Retrieve and plot displacement magnitude
            %   y = obj.Get_Harmonic_Elasticyty_Point('Syy', id_vec); % Retrieve and plot stress in y-direction

            switch field
                case 'u'            
                    % Modulus of the displacement (magnitude)
                    u = obj.femSPACE_elasticity.FIELD_component( obj.Solution.HarmonicElasticity.u , 'norm' );
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Harmonic Displacement $u$ (m)', 'Interpreter', 'latex')
                    y = u;  % Return displacement magnitude for the given region

                case 'ux'            
                    % X-displacement
                    y = obj.femSPACE_elasticity.POSTplotScalar1D( obj.Solution.HarmonicElasticity.ux , id_vec );
                    title('Harmonic Displacement $u_x$ (m)', 'Interpreter', 'latex')

                case 'uy'            
                    % Y-displacement
                    y = obj.femSPACE_elasticity.POSTplotScalar1D( obj.Solution.HarmonicElasticity.uy , id_vec );

                case 'Sxx'    
                    % Stress in the x-direction (sigma_xx)
                    u = obj.Solution.HarmonicElasticity.uSxx;
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{xx}$ (Pa)', 'Interpreter', 'latex') 
                    y = u;  % Return stress in the x-direction for the given region

                case 'Sxy'    
                    % Shear stress (sigma_xy)
                    u = obj.Solution.HarmonicElasticity.uSxy;
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{xy}$ (Pa)', 'Interpreter', 'latex')
                    y = u;  % Return shear stress for the given region

                case 'Syy'    
                    % Stress in the y-direction (sigma_yy)
                    u = obj.Solution.HarmonicElasticity.uSyy;
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec );
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\sigma_{yy}$ (Pa)', 'Interpreter', 'latex')
                    y = u;  % Return stress in the y-direction for the given region
            end
         end

        
        
        function Plot_VibrationModes( obj , field , id_vec , modeNumber )
            % Plot_VibrationModes: This function plots the vibration modes of the system.
            %
            % The function generates plots of the displacement fields corresponding to each of the first
            % `modeNumber` vibration modes for a specified field. It visualizes the mode shapes of the system,
            % including the associated frequency, and outputs a 2D plot for each mode. The field can be any 
            % relevant physical quantity, such as displacement in the x-direction, y-direction, or the magnitude.
            %
            % Inputs:
            %   obj        - The object that contains the solution data for the vibration modes.
            %   field      - A string specifying which field to retrieve and plot. Possible values include:
            %                - 'u'  : Magnitude of displacement
            %                - 'ux' : Displacement in the x-direction
            %                - 'uy' : Displacement in the y-direction
            %   id_vec     - A vector of element IDs specifying the region or mesh area to be plotted.
            %   modeNumber - The number of vibration modes to plot (an integer indicating how many modes to visualize).
            %
            % Outputs:
            %   None. The function generates a series of 2D plots for the selected vibration modes.
            %
            % Process:
            %   - The function loops over the first `modeNumber` vibration modes.
            %   - For each mode, it extracts the corresponding displacement field and plots the 2D shape of the mode.
            %   - The frequency for each mode is printed to the console in Hz.
            %
            % Example Usage:
            %   obj.Plot_VibrationModes('u', id_vec, 3);   % Plot the first 3 vibration modes based on displacement magnitude
            %   obj.Plot_VibrationModes('ux', id_vec, 5);  % Plot the first 5 vibration modes based on displacement in the x-direction

            % Get the current figure number
            figN1 = get(gcf,'Number');

            % Loop through each mode
            for i1 = 1:modeNumber
                % Print the frequency of the current mode
                fprintf( 'Mode f(%d) = %f Hz \n' , i1 , obj.Solution.VibrationModes.vfreqs(i1))

                % Create a new figure for each mode
                figure( figN1 + i1 )

                % Extract the displacement field for the current mode
                us = obj.femSPACE_elasticity.FIELD_component( obj.Solution.VibrationModes.u(:,i1), field );

                % Plot the displacement field for the selected region
                obj.femSPACE_elasticity.POSTplotScalar2D( us(:,1), id_vec )

                % Add labels and title
                xlabel('X(m)');
                ylabel('Y(m)');
                title( sprintf('Mode %d', i1) );
            end
        end

        
        
        function Plot_Estatic_Magnetic( obj , field , id_vec )
            % Plot_Estatic_Magnetic
            %
            % This function generates 2D plots for different static magnetic field quantities.
            % The field to plot is selected based on the 'field' argument. The function supports 
            % various magnetic quantities such as the magnetic vector potential A, the magnetic field B, 
            % the magnetic intensity H, and related components. 
            % 
            % Inputs:
            %   - obj: The object containing the FEM solution, which includes 
            %          the magnetic field solution, relative permeability, and necessary methods.
            %   - field: A string indicating which field to plot. Possible values:
            %     'A'    - Vector potential A (T.m)
            %     'B'    - Magnetic field B (T)
            %     'Bx'   - x-component of magnetic field B (T)
            %     'By'   - y-component of magnetic field B (T)
            %     'Bvec' - Vector of magnetic field B (for GMSH visualization)
            %     'H'    - Magnetic field intensity H (A/m)
            %     'Hx'   - x-component of H (A/m)
            %     'Hy'   - y-component of H (A/m)
            %     'mux'  - Relative permeability in x-direction (dimensionless)
            %     'muy'  - Relative permeability in y-direction (dimensionless)
            %     'Hvec' - Vector of magnetic field intensity H (for GMSH visualization)
            %   - id_vec: Identifier for the region/mesh over which to plot the field (e.g., the finite element mesh).
            %
            % Outputs:
            %   - None. This function creates 2D plots for the selected field and adds labels and titles.

            % Extract the x and y components of the magnetic field (B), magnetization (M),
            % and magnetic intensity (H) from the solution stored in the object 'obj'.
            uBx = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uB , 'x' );  
            uMx = obj.femSPACE_em.FIELD_component( obj.P0FIELDS.uMd , 'x' );
            uHx = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uH , 'x' );
            uBy = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uB , 'y' );  
            uMy = obj.femSPACE_em.FIELD_component( obj.P0FIELDS.uMd , 'y' );
            uHy = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uH , 'y' );
            umux = obj.Solution.StaticMagentic.umux;
            umuy = obj.Solution.StaticMagentic.umuy;
            uH = [uHx; uHy];  % Combine H field components into a vector

            % Select the field to plot based on the 'field' parameter.
            switch field
                case 'A'            
                    % Plot the modulus of the magnetic vector potential A.
                    % This represents the magnetic vector potential in Tesla-meters.
                    u = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uA , 'norm' );
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec ); % Plot the scalar field (A)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic vector potential $A$ (T.m)','Interpreter','latex')

                case 'B'            
                    % Plot the modulus of the magnetic field B.
                    % This represents the total magnetic field in Tesla.
                    u = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uB , 'norm' );
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec ); % Plot the scalar field (B)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic field $B$ (T)','Interpreter','latex')

                case 'Bx'            
                    % Plot the x-component of the magnetic field B.
                    % This represents the magnetic field in the x-direction (Tesla).
                    obj.femSPACE_em.POSTplotScalar2D( uBx , id_vec ); % Plot the scalar field (Bx)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic field $B_x$ (T)','Interpreter','latex')

                case 'By'    
                    % Plot the y-component of the magnetic field B.
                    % This represents the magnetic field in the y-direction (Tesla).
                    obj.femSPACE_em.POSTplotScalar2D( uBy , id_vec ); % Plot the scalar field (By)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic field $B_y$ (T)','Interpreter','latex')  

                case 'Bvec'
                    % Save the magnetic field vector B for visualization in GMSH.
                    % This will generate a vector plot for the magnetic field.
                    obj.femSPACE_em.POSTplotGMSH2Dvector( obj.Solution.StaticMagentic.uB , id_vec , 'Bvec' );

                case 'H'            
                    % Plot the modulus of the magnetic field intensity H.
                    % This represents the magnetic field intensity in A/m.
                    u = obj.femSPACE_em.FIELD_component( obj.Solution.StaticMagentic.uH , 'norm' );
                    obj.femSPACE_elasticity.POSTplotScalar2D( u , id_vec ); % Plot the scalar field (H)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic $H$ field (A/m)','Interpreter','latex')

                case 'Hx'    
                    % Plot the x-component of the magnetic field intensity H.
                    % This represents the magnetic field intensity in the x-direction (A/m).
                    obj.femSPACE_em.POSTplotScalar2D( uHx , id_vec ); % Plot the scalar field (Hx)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic field $H_x$ (A/m)','Interpreter','latex')

                case 'Hy'    
                    % Plot the y-component of the magnetic field intensity H.
                    % This represents the magnetic field intensity in the y-direction (A/m).
                    obj.femSPACE_em.POSTplotScalar2D( uHy , id_vec ); % Plot the scalar field (Hy)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('Static Magnetic field $H_y$ (A/m)','Interpreter','latex')

                case 'mux'    
                    % Plot the relative permeability in the x-direction (dimensionless).
                    % This shows how permeability varies in the x-direction.
                    obj.femSPACE_em.POSTplotScalar2D( umux/obj.mu0 , id_vec ); % Plot the scalar field (mux)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\mu_{x,r}$','Interpreter','latex')

                case 'muy'    
                    % Plot the relative permeability in the y-direction (dimensionless).
                    % This shows how permeability varies in the y-direction.
                    obj.femSPACE_em.POSTplotScalar2D( umuy/obj.mu0 , id_vec ); % Plot the scalar field (muy)
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    title('$\mu_{y,r}$','Interpreter','latex')

                case 'Hvec'
                    % Save the magnetic field intensity vector H for visualization in GMSH.
                    % This will generate a vector plot for the magnetic field intensity.
                    obj.femSPACE_em.POSTplotGMSH2Dvector( uH , id_vec , 'Hvec' );
            end
        end

        
        
        function Plot_DynamicFields( obj , field , id_vec , k , DeformationScale )          
           
            uAz = obj.Solution.Dynamic.uAT(1:end-3,k);
            uBx = obj.femSPACE_em.FIELD_component( obj.Solution.Dynamic.uB(:,k) , 'x' );  
            uBy = obj.femSPACE_em.FIELD_component( obj.Solution.Dynamic.uB(:,k) , 'y' );  
            uBnorm = obj.femSPACE_em.FIELD_component( obj.Solution.Dynamic.uB(:,k) , 'norm' );  
            uMx = obj.femSPACE_em.FIELD_component( obj.P0FIELDS.uMd , 'x' );
            uHx = uBx./obj.Solution.Dynamic.umux(:,k) - uMx;
            uMy = obj.femSPACE_em.FIELD_component( obj.P0FIELDS.uMd , 'y' );
            uHy = uBy./obj.Solution.Dynamic.umuy(:,k) - uMy;
            umux = obj.Solution.Dynamic.umux(:,k);
            umuy = obj.Solution.Dynamic.umuy(:,k);
            uy  = obj.Solution.Dynamic.uy(:,k);

            if ~exist( 'DeformationScale' , 'var' )
                DeformationScale = 1;
            end
            
            %move mesh (magnetic)
                obj.femSPACE_em.msh.POS1(:,1) = obj.femSPACE_em.msh.POS(:,1) + obj.Solution.Dynamic.ux(:,k)*DeformationScale;
                obj.femSPACE_em.msh.POS1(:,2) = obj.femSPACE_em.msh.POS(:,2) + obj.Solution.Dynamic.uy(:,k)*DeformationScale;
                obj.femSPACE_em.prepareMesh(); 
            
            switch field                
                case 'Az'            
                    %modulus of A potential vector
                    obj.femSPACE_em.POSTplotScalar2D( uAz , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic vector potential $A_z$ (T.m) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex')
                case 'B'            
                    %modulus of B magentic field
                    u = obj.femSPACE_em.FIELD_component( obj.Solution.Dynamic.uB(:,k), 'norm' );
                    obj.femSPACE_em.POSTplotScalar2D( u , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic fiel $|B|$ (T) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex')
                case 'Bx'            
                    %Bx                   
                    obj.femSPACE_em.POSTplotScalar2D( uBx , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic fiel $B_x$ (T) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex')
                case 'By'    
                    %By                    
                    obj.femSPACE_em.POSTplotScalar2D( uBy , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic fiel $B_x$ (T) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex')     
                case 'Bnorm'    
                    %By                    
                    obj.femSPACE_em.POSTplotScalar2D( uBnorm , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic fiel $B$ (T) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex') 
                case 'Hx'    
                    %Magnetic field intesnty
                    obj.femSPACE_em.POSTplotScalar2D( uHx , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic fiel $H_x$ (A/m) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex')  
                case 'Hy'    
                    %Magnetic field intesnty
                    obj.femSPACE_em.POSTplotScalar2D( uHy , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title( sprintf( 'Magnetic fiel $H_y$ (A/m) at $t=$ %f s' , obj.tvec(k) ),'Interpreter','latex')  
                    
                case 'mux'    
                    %Magnetic field intesnty
                    obj.femSPACE_em.POSTplotScalar2D( umux/obj.mu0 , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title('$\mu_{x,r}$','Interpreter','latex')
                case 'muy'    
                    %Magnetic field intesnty
                    obj.femSPACE_em.POSTplotScalar2D( umuy/obj.mu0 , id_vec );xlabel('X(m)');ylabel('Y(m)');
                    title('$\mu_{y,r}$','Interpreter','latex')
                case 'uy'
                    %Saves B field vector toplot in gmsh
                    obj.femSPACE_em.POSTplotScalar2D( uy , id_vec );
                    title('$u_y$ displacement (m)','Interpreter','latex')

            end
        end
        
        
       function [y,t] = Plot_Dynamic_Escalars( obj , field )
            % Plot_Dynamic_Escalars
            %
            % This function computes and plots dynamic scalar quantities (voltage, current, power, etc.)
            % from the dynamic solution stored in the 'obj'. It processes the time-dependent data
            % and plots or returns the requested quantity based on the 'field' argument.
            %
            % Inputs:
            %   - obj: The object containing the dynamic FEM solution with time-dependent data.
            %   - field: A string specifying the scalar quantity to plot or compute. Possible values:
            %       'V'     - Total voltage (sum of Vin and Vout) over time.
            %       'Vin'   - Input voltage over time.
            %       'Vout'  - Output voltage over time.
            %       'I'     - Current over time.
            %       'Vrms'  - Root mean square (RMS) value of the total voltage.
            %       'Irms'  - RMS value of the current.
            %       'Peff'  - Effective power calculated using RMS current and resistance RL.
            %       'Pdamp' - Damping loss power calculated over all domains.
            %
            % Outputs:
            %   - y: The computed value or vector of values corresponding to the requested field.
            %   - t: The corresponding time vector for the plotted or returned data.

            % Preallocate arrays for the computed quantities
            Vin = zeros(1, size(obj.Solution.Dynamic.uAT, 2));
            Vout = zeros(1, size(obj.Solution.Dynamic.uAT, 2));
            I = zeros(1, size(obj.Solution.Dynamic.uAT, 2));
            t = zeros(1, size(obj.Solution.Dynamic.uAT, 2));

            % Compute Vin, Vout, I, and corresponding time vector t by averaging adjacent time steps
            for k = 2:size(obj.Solution.Dynamic.uAT, 2)
                Vin(k)  =  ( obj.Solution.Dynamic.uAT(end-2, k) + obj.Solution.Dynamic.uAT(end-2, k-1) ) / 2;
                Vout(k) =  ( obj.Solution.Dynamic.uAT(end-1, k) + obj.Solution.Dynamic.uAT(end-1, k-1) ) / 2;
                I(k)    =  ( obj.Solution.Dynamic.uAT(end, k) + obj.Solution.Dynamic.uAT(end, k-1) ) / 2;
                t(k)    =  ( obj.tvec(k) + obj.tvec(k-1) ) / 2;
            end

            % Initialize the first elements as zero
            Vin(1) = 0;
            Vout(1) = 0;
            I(1) = 0;

            % Compute total voltage V as the sum of Vin and Vout
            V = Vin + Vout;

            % Switch case to handle different field types
            switch field    
                case 'V'                    
                    % Plot total voltage over time
                    plot(t, V); 
                    xlabel('$t$ (s)', 'Interpreter', 'latex'); 
                    ylabel('$V$', 'Interpreter', 'latex')
                    y = V;

                case 'Vin'                    
                    % Plot input voltage over time
                    plot(t, Vin); 
                    xlabel('$t$ (s)', 'Interpreter', 'latex'); 
                    ylabel('$V_{in}$', 'Interpreter', 'latex')
                    y = Vin;

                case 'Vout'                    
                    % Plot output voltage over time
                    plot(t, Vout); 
                    xlabel('$t$ (s)', 'Interpreter', 'latex'); 
                    ylabel('$V_{out}$', 'Interpreter', 'latex')
                    y = Vout;

                case 'I'                    
                    % Plot current over time
                    plot(t, I); 
                    xlabel('$t$ (s)', 'Interpreter', 'latex'); 
                    ylabel('$I$ (A)', 'Interpreter', 'latex')
                    y = I;

                case 'Vrms'                    
                    % Compute RMS value of voltage over the last 'obj.Steps' time steps
                    y = mean(sqrt(V(end - obj.Steps:end).^2));

                case 'Irms'                    
                    % Compute RMS value of current over the last 'obj.Steps' time steps
                    y = mean(sqrt(I(end - obj.Steps:end).^2));

                case 'Peff'                    
                    % Compute effective power using RMS current and resistance RL
                    Irms = mean(sqrt(I(end - obj.Steps:end).^2));
                    y = obj.RL * Irms^2;

                case 'Pdamp'
                    % Initialize power damping
                    y = 0;

                    % Loop over each unique domain to compute damping loss power
                    for iD = unique(obj.femSPACE_em.msh.TRIANGLES(:,end))'
                        % Project P0 field to P1 for rho_alpha calculation
                        uP1_rho_alpha = obj.femSPACE_em.POSTprojectP0toP1( ...
                            obj.alpha * obj.P0FIELDS.urho, iD, 'SurfaceIntegral');

                        % Compute total displacement
                        uP1_u2 = obj.Solution.HarmonicElasticity.ux .* conj(obj.Solution.HarmonicElasticity.ux) + ...
                                 obj.Solution.HarmonicElasticity.uy .* conj(obj.Solution.HarmonicElasticity.uy);

                        % Damping loss power density calculation
                        uP1_rhoD = 1/2 * obj.omega^2 * uP1_rho_alpha .* uP1_u2;

                        % Compute integral over the surface
                        [ ~ , ysum ] = obj.femSPACE_em.POSTmean( 2, iD, 0.5 * obj.omega^2 * (uP1_rho_alpha .* uP1_u2));

                        % Add damping loss power from this domain to the total
                        y = y + ysum;
                    end

                    % Convert surface integral to volume by multiplying with object thickness
                    y = y * obj.Lz;
            end       
        end

        
        
        function Plot_Dynamic( obj , field , id )
            % Plot_Dynamic
            %
            % This function processes and plots dynamic quantities (magnetic, mechanical, etc.)
            % from the simulation results stored in the 'obj'. It computes the average values of
            % various quantities over time for a given region specified by 'id'.
            %
            % Inputs:
            %   - obj: The object containing the dynamic FEM solution data.
            %   - field: A string specifying the quantity to plot. Possible values include:
            %       'muxr'   - Relative permeability in the x-direction.
            %       'Hx'     - Magnetic field intensity in the x-direction (A/m).
            %       'Bx'     - Magnetic flux density in the x-direction (T).
            %       'Bnorm'  - Norm of the magnetic flux density (T).
            %       'Sxx'    - Mechanical stress in the x-direction (Pa).
            %       'ux'     - Displacement in the x-direction (m).
            %       'uy'     - Displacement in the y-direction (m).
            %       'Fmagx'  - Magnetic force component in the x-direction (N).
            %       'Fmagy'  - Magnetic force component in the y-direction (N).
            %   - id: The region identifier used for computing average values over surfaces.
            %
            % Outputs:
            %   - A plot of the requested field over time.

            % Preallocate arrays for storing computed values
            numSteps = length(obj.tvec);
            t = zeros(1, numSteps);
            mux = zeros(1, numSteps);
            muy = zeros(1, numSteps);
            Hx = zeros(1, numSteps);
            Bx = zeros(1, numSteps);
            Bnorm = zeros(1, numSteps);
            Sxx = zeros(1, numSteps);
            ux = zeros(1, numSteps);
            uy = zeros(1, numSteps);
            uFmagx = zeros(1, numSteps);
            uFmagy = zeros(1, numSteps);

            % Loop over all time steps to compute mean values for each field
            for k = 1:numSteps
                t(k)       = obj.tvec(k);
                mux(k)     = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.Solution.Dynamic.umux(:,k)); 
                muy(k)     = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.Solution.Dynamic.umuy(:,k));
                Hx(k)      = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.femSPACE_em.FIELD_component(obj.Solution.Dynamic.uH(:,k), 'x'));
                Bx(k)      = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.femSPACE_em.FIELD_component(obj.Solution.Dynamic.uB(:,k), 'x'));
                Bnorm(k)   = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.femSPACE_em.FIELD_component(obj.Solution.Dynamic.uB(:,k), 'norm'));
                Sxx(k)     = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.Solution.Dynamic.uSxx(:,k));
                ux(k)      = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.Solution.Dynamic.ux(:,k));
                uy(k)      = obj.femSPACE_em.POSTmean('SurfaceIntegral', id, obj.Solution.Dynamic.uy(:,k));
            end

            % Plotting the requested field
            switch field
                case 'muxr'                    
                    plot(t, mux/obj.mu0);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<\mu_{r,x}>$', 'Interpreter', 'latex');

                case 'Hx'
                    plot(t, Hx);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<H_x>$ (A/m)', 'Interpreter', 'latex');

                case 'Bx'
                    plot(t, Bx);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<B_x>$ (T)', 'Interpreter', 'latex');

                case 'Bnorm'
                    plot(t, Bnorm);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<B>$ (T)', 'Interpreter', 'latex');

                case 'Sxx'
                    plot(t, Sxx);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<S_{xx}>$ (Pa)', 'Interpreter', 'latex');

                case 'ux'
                    plot(t, ux);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<u_{x}>$ (m)', 'Interpreter', 'latex');

                case 'uy'
                    plot(t, uy);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$<u_{y}>$ (m)', 'Interpreter', 'latex');

                case 'Fmagx'
                    plot(t, uFmagx);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$F_{mag,x}$ (N)', 'Interpreter', 'latex');

                case 'Fmagy'
                    plot(t, uFmagy);
                    xlabel('$t$ (s)', 'Interpreter', 'latex');
                    ylabel('$F_{mag,y}$ (N)', 'Interpreter', 'latex');
            end       
        end

        
        
        
        function Plot_DynamicFieldsGMSH( obj , field , id_vec , k , DeformationScale )
            % Plot_DynamicFieldsGMSH
            %
            % This function generates plots or saves field data related to magnetic
            % vector potentials, magnetic fields, and forces using GMSH-compatible formats.
            % The mesh can be deformed according to the dynamic displacement fields.
            %
            % Inputs:
            %   - obj: The object containing the FEM simulation data.
            %   - field: A string specifying the field to be plotted or saved. Possible values include:
            %       'A'      - Magnetic vector potential magnitude.
            %       'Fvec'   - Magnetic force field vector (saved for GMSH visualization).
            %       'Bx'     - Magnetic field component Bx (T).
            %       'By'     - Magnetic field component By (T).
            %       'Bvec'   - Magnetic field vector (saved for GMSH visualization).
            %       'Hx'     - Magnetic field intensity component Hx (A/m).
            %       'Hy'     - Magnetic field intensity component Hy (A/m).
            %       'mux'    - Relative permeability in the x-direction.
            %       'muy'    - Relative permeability in the y-direction.
            %       'Hvec'   - Magnetic field vector (H) (saved for GMSH visualization).
            %   - id_vec: Indices of the regions to be plotted.
            %   - k: The time step index.
            %   - DeformationScale: Optional. A scaling factor for mesh deformation (default: 1).
            %
            % Outputs:
            %   - Displays or saves plots depending on the selected field.

            % Extract field components
            uBx = obj.femSPACE_em.FIELD_component(obj.Solution.Dynamic.uB(:,k), 'x');  
            uMx = obj.femSPACE_em.FIELD_component(obj.P0FIELDS.uMd, 'x');            
            uBy = obj.femSPACE_em.FIELD_component(obj.Solution.Dynamic.uB(:,k), 'y');  
            uMy = obj.femSPACE_em.FIELD_component(obj.P0FIELDS.uMd, 'y');
            uHx = uBx ./ obj.Solution.Dynamic.umux(:,k) - uMx;
            uHy = uBy ./ obj.Solution.Dynamic.umuy(:,k) - uMy;
            umux = obj.Solution.Dynamic.umux(:,k);
            umuy = obj.Solution.Dynamic.umuy(:,k);
            uH = [uHx; uHy];

            % Set default deformation scale if not provided
            if ~exist('DeformationScale', 'var')
                DeformationScale = 1;
            end

            % Deform the mesh according to the displacement fields (if applicable)
            obj.femSPACE_em.msh.POS1(:,1) = obj.femSPACE_em.msh.POS(:,1) + obj.Solution.Dynamic.ux(:,k) * DeformationScale;
            obj.femSPACE_em.msh.POS1(:,2) = obj.femSPACE_em.msh.POS(:,2) + obj.Solution.Dynamic.uy(:,k) * DeformationScale;
            obj.femSPACE_em.prepareMesh(); 

            % Handle different cases for the requested field
            switch field                
                case 'A'
                    % Plot magnitude of magnetic vector potential A
                    u = obj.femSPACE_em.FIELD_component(obj.Solution.Dynamic.uAT(1:end-3, k), 'norm');
                    obj.femSPACE_elasticity.POSTplotScalar2D(u, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title(sprintf('Magnetic vector potential $A$ (T.m) at $t=$ %f s', obj.tvec(k)), 'Interpreter', 'latex');

                case 'Fvec'
                    % Save magnetic force vector for GMSH visualization
                    obj.femSPACE_em.POSTplotGMSH2Dvector(obj.Solution.Dynamic.uFmag(:,k), id_vec, sprintf('Fvec_t%fs', obj.tvec(k)));

                case 'Bx'
                    % Plot Bx field
                    obj.femSPACE_em.POSTplotScalar2D(uBx, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title(sprintf('Magnetic field $B_x$ (T) at $t=$ %f s', obj.tvec(k)), 'Interpreter', 'latex');

                case 'By'
                    % Plot By field
                    obj.femSPACE_em.POSTplotScalar2D(uBy, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title(sprintf('Magnetic field $B_y$ (T) at $t=$ %f s', obj.tvec(k)), 'Interpreter', 'latex');  

                case 'Bvec'
                    % Save magnetic field vector for GMSH visualization
                    obj.femSPACE_em.POSTplotGMSH2Dvector(obj.Solution.Dynamic.uB(:,k), id_vec, sprintf('Bvec_t%fs', obj.tvec(k)));

                case 'Hx'
                    % Plot Hx field
                    obj.femSPACE_em.POSTplotScalar2D(uHx, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title(sprintf('Magnetic field $H_x$ (A/m) at $t=$ %f s', obj.tvec(k)), 'Interpreter', 'latex');

                case 'Hy'
                    % Plot Hy field
                    obj.femSPACE_em.POSTplotScalar2D(uHy, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title(sprintf('Magnetic field $H_y$ (A/m) at $t=$ %f s', obj.tvec(k)), 'Interpreter', 'latex');

                case 'mux'
                    % Plot relative permeability in the x-direction
                    obj.femSPACE_em.POSTplotScalar2D(umux / obj.mu0, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title('$\mu_{x,r}$', 'Interpreter', 'latex');

                case 'muy'
                    % Plot relative permeability in the y-direction
                    obj.femSPACE_em.POSTplotScalar2D(umuy / obj.mu0, id_vec);
                    xlabel('X(m)'); ylabel('Y(m)');
                    title('$\mu_{y,r}$', 'Interpreter', 'latex');

                case 'Hvec'
                    % Save magnetic field vector (H) for GMSH visualization
                    obj.femSPACE_em.POSTplotGMSH2Dvector(uH, id_vec, 'Hvec');
            end
        end

        
        
                
        
        function Plot_Magnetostrictive( obj , number , Hmax , Smax )
            % Plot_Magnetostrictive
            %
            % This function plots the magnetization (M) and relative permeability (mu_r)
            % of a specified magnetostrictive material from the object's property data.
            % Two 3D surface plots are generated:
            %   - Left: Magnetization M as a function of Magnetic Field H and Stress Sxx.
            %   - Right: Relative Permeability (mu_r) as a function of H and Sxx.
            %
            % Inputs:
            %   - obj: The object containing the magnetostrictive properties.
            %   - number: Index of the material to be plotted from obj.Prop_Magnetostrictive.
            %   - Hmax (optional): Maximum value of H to be displayed (default: obj.Hmax).
            %   - Smax (optional): Maximum value of Sxx to be displayed (default: obj.Smax).
            %
            % Outputs:
            %   - Two subplots illustrating the Magnetization and Relative Permeability.

            % Set default values for Hmax and Smax if not provided
            if nargin < 3, Hmax = obj.Hmax; end
            if nargin < 4, Smax = obj.Smax; end

            % Extract material properties for the given number index
            HH = obj.Prop_Magnetostrictive{number}.HH;
            SS = obj.Prop_Magnetostrictive{number}.SS;
            MM = obj.Prop_Magnetostrictive{number}.MM;
            MUMU = obj.Prop_Magnetostrictive{number}.MUMU;

            % Plot 1: Magnetization (M) as a function of H and Sxx
            subplot(1,2,1);
            surf(HH, SS, MM);
            xlabel('H (A/m)');
            ylabel('$\sigma_{xx}$ (Pa)', 'Interpreter', 'latex');
            zlabel('$\mu_0 M$ (T)', 'Interpreter', 'latex');
            title('Magnetization', 'Interpreter', 'latex');
            axis([0 Hmax -Smax Smax]);     

            % Plot 2: Relative Permeability (mu_r) as a function of H and Sxx
            subplot(1,2,2);
            surf(HH, SS, MUMU / obj.mu0);
            xlabel('H (A/m)');
            ylabel('$\sigma_{xx}$ (Pa)', 'Interpreter', 'latex');
            zlabel('$\mu_r$', 'Interpreter', 'latex');
            title('Permeability $\mu_r$', 'Interpreter', 'latex');
            axis([0 Hmax -Smax Smax]);
        end
    
        
        
        
        function Plot_NoLinearPermeability( obj , number , Hmax )
            % Plot_NoLinearPermeability - Plots the magnetization curve and 
            % relative permeability curve for a nonlinear permeability material.
            %
            % Syntax:
            %   obj.Plot_NoLinearPermeability(number, Hmax)
            %
            % Inputs:
            %   obj    - The object containing the nonlinear permeability properties.
            %   number - Index of the nonlinear permeability property to plot.
            %   Hmax   - Maximum value of applied magnetic field H for x-axis limit.
            %
            % Description:
            %   This function generates two subplots:
            %     1. Magnetization Curve: Plots the applied magnetic field (HH) against 
            %        the magnetic flux density (BB).
            %     2. Relative Permeability Curve: Plots the applied magnetic field (HH) 
            %        against the relative permeability (MUMU / obj.mu0).
            %
            % Example:
            %   obj.Plot_NoLinearPermeability(1, 1000);
            %
            % See also:
            %   None

            subplot(1,2,1);title('Magnetization')
            HH = obj.Prop_NoLinearPermeability{number}.HH;               
            BB = obj.Prop_NoLinearPermeability{number}.BB;
            MUMU = obj.Prop_NoLinearPermeability{number}.MUMU;
            
            % Plot Magnetization Curve
            plot( HH , BB );
            xlabel('$H_{app}$(A/m)','Interpreter','latex');
            ylabel('$B$(T)','Interpreter','latex');
            xlim([-Hmax, Hmax]);
            
            subplot(1,2,2);title('permeavility $\mu_r$','Interpreter','latex')
            
            % Plot Relative Permeability Curve
            plot( HH , MUMU/obj.mu0 );
            xlabel('H(A/m)');
            ylabel('$\mu_r$','Interpreter','latex')
            xlim([-Hmax, Hmax]);
        end

        
        
        function Plot_Surface_RegionID( obj )
            % Plot_Surface_RegionID - Plots the surface region IDs of a mesh.
            %
            % Syntax:
            %   obj.Plot_Surface_RegionID()
            %
            % Inputs:
            %   obj - The object containing the mesh and plotting functionalities.
            %
            % Description:
            %   This function plots a 2D scalar field representing region IDs of the 
            %   surface mesh elements. The IDs are extracted from the 'TRIANGLES' 
            %   matrix and plotted using the POSTplotScalar2D method.
            %
            % Example:
            %   obj.Plot_Surface_RegionID();
            %
            % See also:
            %   POSTplotScalar2D

            u = obj.femSPACE_elasticity.msh.TRIANGLES(:,4);
            obj.femSPACE_elasticity.POSTplotScalar2D( u , 0 );
            xlabel('X(m)');
            ylabel('Y(m)');
        end      
       
    end
end

