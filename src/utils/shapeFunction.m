% Copyright (C) 2005 Garth N .Wells
%
% Compute shape functions and determinant of the Jacobian
%
%

% Format
%
% N(i)    = shape functions of node i
% dN(i,j) = derivative of shape functions i in the direction x_j
% j       = determinant of the Jacobian
%

function [N, dN, j] = shapeFunction(x, gp, type)
  

  switch(type) 
    case 'Bar1'
      N = zeros(2,1); dN = zeros(2,1); J = zeros(1,1); 

      xi  = gp(2);

      N(1) = -xi/2 + 0.5;
      N(2) =  xi/2 + 0.5;

      dN(1,1) = -1/2; 
      dN(2,1) =  1/2;
 
    case 'Bar2'

      N = zeros(3,1); dN = zeros(3,1); J = zeros(1,1); 

      xi  = gp(2);

      N(1) =  xi*xi/2 - xi/2;
      N(2) =  xi*xi/2 + xi/2;
      N(3) = -xi*xi + 1;

      dN(1,1) =  xi - 0.5;
      dN(2,1) =  xi + 0.5;
      dN(3,1) = -2*xi;
 
	  case 'Tri3' 

      N = zeros(3,1); dN = zeros(3,2); J = zeros(2,2);

      xi  = gp(2);
      eta = gp(3);

      N(1)   = 1.0 - xi - eta;
      N(2)   = xi;
      N(3)   = eta;

      dN(1,1) = -1.0;
      dN(1,2) = -1.0;
      dN(2,1) =  1.0;
      dN(2,2) =  0.0;
      dN(3,1) =  0.0;
      dN(3,2) =  1.0;

	  case 'Tri6' 

      disp('Tri6 not yet programmed ')
 
    case 'Quad4'  % shape functions for 4-noded quad

		  N = zeros(4,1); dN = zeros(4,2); J = zeros(2,2); 
				
		  a(1) =-1.0; a(2) = 1.0; a(3) = 1.0; a(4) =-1.0; % x locations of nodes
		  b(1) =-1.0; b(2) =-1.0; b(3) = 1.0; b(4) = 1.0; % y locations of nodes
		
		  xi  = gp(2);
		  eta = gp(3);
		
		  N(:)    = 0.25*(1.0 + a(:)*xi + b(:)*eta + a(:).*b(:)*xi*eta);

  		dN(:,1) = 0.25*(a(:) + a(:).*b(:)*eta);
	   	dN(:,2) = 0.25*(b(:) + a(:).*b(:)*xi);

	  case 'Quad8'

      N = zeros(8,1); dN = zeros(8,2); J = zeros(2,2); 

      a(1) =-1.0; a(2) = 1.0; a(3) = 1.0; a(4) =-1.0; % x locations of nodes
      b(1) =-1.0; b(2) =-1.0; b(3) = 1.0; b(4) = 1.0; % y locations of nodes

      xi  = gp(2);
      eta = gp(3);

      % Compute shape functions for Quad4 element
      N(1:4) = 0.25*(1.0 + a(:)*xi + b(:)*eta + a(:).*b(:)*xi*eta);

      % Compute mid-side node shape functions for Quad8
      N(5)   = 0.5*(1.0-xi^2)*(1.0-eta);
      N(6)   = 0.5*(1.0-eta^2)*(1.0+xi);
      N(7)   = 0.5*(1.0-xi^2)*(1.0+eta);
      N(8)   = 0.5*(1.0-eta^2)*(1.0-xi);

      % Modify corner nodes to make quadratic
      N(1)   = N(1) - 0.5*(N(5) + N(8));
      N(2)   = N(2) - 0.5*(N(5) + N(6));
      N(3)   = N(3) - 0.5*(N(6) + N(7));
      N(4)   = N(4) - 0.5*(N(7) + N(8));

      % Compute derivatives
      dN(1:4,1) = 0.25*(a(:) + a(:).*b(:)*eta);
      dN(1:4,2) = 0.25*(b(:) + a(:).*b(:)*xi);

      dN(5,1)   = 0.5*(-2.0*xi)*(1.0-eta);
      dN(5,2)   = 0.5*(1.0-xi^2)*(-1.0);
      dN(6,1)   = 0.5*(1.0-eta^2)*(1.0);
      dN(6,2)   = 0.5*(-2.0*eta)*(1.0+xi);
      dN(7,1)   = 0.5*(-2.0*xi)*(1.0+eta);
      dN(7,2)   = 0.5*(1.0-xi^2)*(1.0);
      dN(8,1)   = 0.5*(1.0-eta^2)*(-1.0);
      dN(8,2)   = 0.5*(-2.0*eta)*(1.0-xi);

      dN(1,1) = dN(1,1) - 0.5*(dN(5,1) + dN(8,1));
      dN(1,2) = dN(1,2) - 0.5*(dN(5,2) + dN(8,2));
      dN(2,1) = dN(2,1) - 0.5*(dN(5,1) + dN(6,1));
      dN(2,2) = dN(2,2) - 0.5*(dN(5,2) + dN(6,2));
      dN(3,1) = dN(3,1) - 0.5*(dN(6,1) + dN(7,1));
      dN(3,2) = dN(3,2) - 0.5*(dN(6,2) + dN(7,2));
      dN(4,1) = dN(4,1) - 0.5*(dN(7,1) + dN(8,1));
      dN(4,2) = dN(4,2) - 0.5*(dN(7,2) + dN(8,2));

	  case 'Quad9' 

        disp('Quad9 not yet programmed ')

	  case 'Tet4' 

        disp('Tet4 not yet programmed ')

	  otherwise 

        disp('Element type not yet programmed ')

	end

  %	transform derivatives from isoparametric configuration to global config.
  
	% form Jacobian matrix
	J = x * dN;  % Note that J is the inverse of the usual J

	% determinant of Jacobain matrix
	j = det(J);
	
  if(j < 0) 
    disp('Warning: negative Jacobian')
  end

	% transform derivatives to global system
	% dN = dN * inv(J);

    dN = dN / J;

    N  = N';
    dN = dN';
end