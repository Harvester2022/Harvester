clear all
close all

%Clase prmt con los parametros wm1 y wm2. Tamaino imanes cuadrados para
%pasar al modelo
prmt.wm1 = 4;
prmt.wm2 = 4;


% Cargar problema y definir geometria

%Creamos el problema FEM harvester
escala = 1e-3; %multiplicar por 1e-3 para pasar de mm a metros. Gmsh esta en mm
gmsh_file = 'Cantilebel_1.geo'; %fichero con la geometria del harvester gmsh para cargar. Abrir fichero con gmsh para ver dominios. 
Harvester1 = MagnetostrictiveHarvester2D( gmsh_file , escala , prmt );

%Propiedades de la simulacion
    %Anchura del problema (perpendicular al plano OXY)
    Harvester1.Lz = 5e-3;
    Harvester1.Periods = 3;
    Harvester1.Steps = 50;    
    
    Harvester1.id_infty = 1;
    
%Propiedades vibraciones aplicadas
    freq = 60.1; %frecuencia en Hz
    az = 1; %amplitud aceleracion en g (x 9.81 m/s^2)
    Harvester1.az = az;    
    Harvester1.omega = freq*2*pi; %frecuencia angular
    
    

% Cargar las porpiedades de los materiales, imanes,...

    %Factores volumen, para multiplicar la imanación volumetrica de los
    %imanes
        xm1 = 5e-3; %Anchura real iman 1 (perpendicular al plano OXY, en metros)
        xm2 = 5e-3; %Anchura real iman 2 (perpendicular al plano OXY, en metros)
        %Imanación de cada iman (A/m) Ferrita?
        Md1 = 400e3;
    
        fV1 = xm1/Harvester1.Lz;
        fV2 = xm2/Harvester1.Lz;
        %Imanacion dos imanes corregida por el factor volumen
        Md1 = Md1*fV1;
        Md2 = Md1*fV2;
        
    %Propiedades mecanicas
        %propiedad elasticidad numero 1: imanes (3,4), Fe (2) frame y GaFe (5).
        Harvester1.Prop_Elasticity{1}.id = [2,3,4,5];
        Harvester1.Prop_Elasticity{1}.YoungModulus = 90e9; %Pa
        Harvester1.Prop_Elasticity{1}.PoissonCoefficient = 0.2;
        Harvester1.Prop_Elasticity{1}.Density = 7874; %(kg/m^3)
        %propiedad elasticidad numero 2: bobina. Material elastico blando
        Harvester1.Prop_Elasticity{2}.id = [6,7];
        Harvester1.Prop_Elasticity{2}.YoungModulus = 50e9;%Pa
        Harvester1.Prop_Elasticity{2}.PoissonCoefficient = 0.2;
        Harvester1.Prop_Elasticity{2}.Density = 500;%(kg/m^3)
        
        %iD de la linea fijada al shaker (u=0)
        Harvester1.id_fix = 2;
        
        %Damping factor alfa, igual en todo el sistema
        Harvester1.alpha = 1; %s^-1
        
        %Mass tip
        m_ferrita = 4e-3;%masa iman de ferrita
        Harvester1.MassTips{1}.id = 3;
        Harvester1.MassTips{1}.value = m_ferrita;
        Harvester1.MassTips{2}.id = 4;
        Harvester1.MassTips{2}.value = m_ferrita;
        
    %Propiedades electricas
        %Conductividad electrica frame de Fe
        Harvester1.Prop_Econdictivity{1}.id = 2;
        Harvester1.Prop_Econdictivity{1}.value = 10e6; %S/m
        %Porpiedades bobinas captadoras
        Harvester1.Coil.In.id = 6;
        Harvester1.Coil.Out.id = 7;
        Harvester1.Coil.Nturn  = 3500;
        Harvester1.Coil.WireSection = pi*(18e-6)^2; % radio de 18 micras
        Harvester1.RL = 1000e3;  %Resistencia de carga RL (ohmios)
        
        %iD de la linea donde el campo magnetico es nulo (Az->0)
        Harvester1.id_infty = 1;

    %Propiedades Magneticas
        %Iman 1
        Harvester1.Magnets{1}.id = 3;
        Harvester1.Magnets{1}.value = [0,+Md1];
        %Iman 2
        Harvester1.Magnets{2}.id = 4;
        Harvester1.Magnets{2}.value = [0,-Md2];
        %Permeabilidad Fe cosntante linear
        Harvester1.Prop_Permeability{1}.id = 2;
        Harvester1.Prop_Permeability{1}.murx = 500;
        Harvester1.Prop_Permeability{1}.mury = 500;
        %Permeabilidad Fe no linear. Sobreescribe la permeavilidad
        %constnate linear definida anteriormente(500)
        Harvester1.Prop_NoLinearPermeability{1}.id = 2;
        Harvester1.Prop_NoLinearPermeability{1}.datafile = 'FeFrame.txt';
        %GaFe lamina
        Harvester1.Prop_Magnetostrictive{1}.id = 5;
        Harvester1.Prop_Magnetostrictive{1}.datafile = 'Galfenol.txt';
        

        
% Solver
    %Cargamos el problema con todo lo definido anteriormente.
    Harvester1.LoadHarvester2DProblem();
    
    %Resolver la elasticidad estatica (efecto gravedad)
    Harvester1.Solve_Estatic_Elasticy();
    
    %Resolver la elasticidad armonica frecuencia freq
    Harvester1.Solve_Harmonic_Elasticyty( );
    
    %Resolver los modos de vibracion.
    Harvester1.Solve_VibrationModes();

    %Reolver la imanación estatica 
    Harvester1.Solve_StaticMagnetic_noLinear( );
    

    
% POST. Plotear. 
    %Escribir las frecuencias de vibracion
    for i=1:length( Harvester1.Solution.VibrationModes.vfreqs )
        fprintf( 'frecuencia Resonancia %d a %f Hz.\n' , i , Harvester1.Solution.VibrationModes.vfreqs(i) ) 
    end
    
    %Plotear las propiedades del GaFe cargadas
        figure(1)
        Harvester1.Plot_Magnetostrictive(1,1e4,2e8);     
        
        figure(2)
        Harvester1.Plot_NoLinearPermeability( 1 , 1e4 )

    %plotear campos 2D
        %Desplacamiento vertical
        %Estatico
        figure(10)
        Harvester1.Plot_Estatic_Elasticyty( 'uy' , [2,3,4,5,6,7] )        
        %Dinamico
        figure(11)
        Harvester1.Plot_Harmonic_Elasticyty( 'uy' , [2,3,4,5,6,7] )
        %tensiones en x e y (no cicalladura Sxy) en el GaFe
        figure(20)
        subplot(1,2,1)
        Harvester1.Plot_Harmonic_Elasticyty( 'Sxx' , [5] )
        subplot(1,2,2)
        Harvester1.Plot_Harmonic_Elasticyty( 'Syy' , [5] )
        %Plotear campo magnetico B estatico (inicial). En el Fe, GaFe e imanes (2,3,4,5)
        figure(30)
        Harvester1.Plot_Estatic_Magnetic( 'B' , [2,3,4,5] )
        %Plotear campo magnetico H estatico (inicial). En el Fe, GaFe e imanes (2,3,4,5)
        figure(30)
        Harvester1.Plot_Estatic_Magnetic( 'H' , [5] )
        
        

        
        
%% Integracion temporal
    %flog: mover mallado. Al no ser electromagnetico (U-shape), todo estatico
    flag_move_mesh = 0;
    Harvester1.Solve_Dynamic( flag_move_mesh );
        
% POST resultados de la integracion tenporal            
    %Plotear evoluciones temporales Solo si hemos integrado tenporalemente, si
    %no, da error, pues no existe el campo 'V' en funcion del tiempo  
        %Plotear el V frente a t. 
        figure(40)
        Harvester1.Plot_Dynamic_Escalars( 'V' );
        
        %Plotear promedio de muxr frente a t en el GaFe (5)
        figure(41)
        Harvester1.Plot_Dynamic( 'muxr' , 5 ) 
        
        %Plotear promedio de Bx (T) frente a t en el GaFe (5)
        figure(42)
        Harvester1.Plot_Dynamic( 'Bx' , 5 )
    
        %Plotear promedio de Sxx (tension) frente a t en el GaFe (5)
        figure(43)
        Harvester1.Plot_Dynamic( 'Sxx' , 5 ) 
        
        %imprimir cual es la potencia y el voltaje effectivos 
        Peff = Harvester1.Plot_Dynamic_Escalars( 'Peff' );
        Vrms = Harvester1.Plot_Dynamic_Escalars( 'Vrms' );
        fprintf( 'Voltaje rms = %e y potencia efectiva Peff = %e  uw\n' , Vrms , Peff*1e6 )
  

    