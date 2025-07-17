function Cijkl = ENSEMBLE_strain_C_tensor( dim , Y , P )
%ENSEMBLE_strain_C_tensor: Ensembles the linear elasticity Cijkl tensor. 
%INPUT
%	- dim: field dimensions (2 or 3). 
%   - Y: scalar P0 field with the Young modulus value in each element. If a single number is provided, a constant value for all domain is assumed. 
%   - P: scalar P0 field with the Poisson coefficient value in each element. If a single number is provided, a constant value for all domain is assumed. 
%OUTPUT: 
	% 	- Cijkl: dim x dim x dim x dim cell.



   Cijkl = cell(dim , dim , dim , dim );
   
   %lamee parameters
   lambda = Y.*P./( 1 + P )./(1-2*P);
   G = Y./2./(1+P);
   
   for i = 1:dim
       for j=1:dim
           for k=1:dim
               for l = 1:dim
                   
                   Cijkl{i,j,k,l}= lambda*(i==k)*(j==l) + G*( (i==j)*(k==l) + (i==l)*(k==j) );
                   
               end
               
           end
           
       end
       
   end