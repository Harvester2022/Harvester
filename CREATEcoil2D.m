function coil = CREATEcoil2D( obj , coil )
%ENSEMBLEcoil2D See documentation

    %FEM space
    femSPACE2D = obj.femSPACE_em;
    
    %If  coil.Out is not defined
    if ~isfield(coil, 'Out')
        coil.Out.id = coil.In.id;    
        coil.Out.wireSection = coil.In.wireSection;
    end
    
    %itegral operators: 
        INTin = femSPACE2D.createVector( 'SurfaceIntegral' , coil.In.id )';
        INTout = femSPACE2D.createVector( 'SurfaceIntegral' , coil.Out.id )';
        
        %Terminal sections 
        coil.In.S  = full( sum( INTin ) );            
        coil.Out.S = full( sum( INTout ) );
        
        %wire densities
        jtin  = - coil.Nturn /coil.In.S;
        jtout = + coil.Nturn /coil.Out.S;

        %Effective average conductivities
        coil.In.sigmaEFF  = coil.WireSection*coil.Nturn*coil.sigma/coil.In.S;
        coil.Out.sigmaEFF = coil.WireSection*coil.Nturn*coil.sigma/coil.Out.S;
        
        %CoilResistances
        coil.In.Rdc = coil.Nturn*obj.Lz/( coil.WireSection*coil.sigma );   
        coil.Out.Rdc = coil.Nturn*obj.Lz/( coil.WireSection*coil.sigma ); 

        %Save in coil class
        JsIN = -jtin/coil.In.Rdc;
        JsOUT= -jtout/coil.Out.Rdc;
        coil.In.BJs  = femSPACE2D.createVector( 'SurfaceIntegral' , coil.In.id  , JsIN  );
        coil.Out.BJs = femSPACE2D.createVector( 'SurfaceIntegral' , coil.Out.id , JsOUT );

        coil.In.emfIntegral    = +INTin*jtin*obj.Lz;
        coil.Out.emfIntegral   = +INTout*jtout*obj.Lz;
        
        coil.uJ = obj.femSPACE_em.FIELD_constantP0( 'Surface' , coil.In.id , jtin  ) + obj.femSPACE_em.FIELD_constantP0( 'Surface' , coil.Out.id , jtout  ) ;

