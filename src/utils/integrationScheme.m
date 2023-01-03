% Copyright (C) 2005 Garth N .Wells
%
% Return integration point locations and weights 
%
% Format
%
% gp(1,j)   = integration weight for point j
% gp(2,j)   = \xi location for point j
% gp(3,j)   = \eta location for point j (2D & 3D)
% gp(4,j)   = \zi location for point j (3D only)

function gp = integrationScheme(points, dim)

  if dim == 1

    if points == 1	% 1-point Gauss integration in 1D

      gp     = zeros(2,1);

      gp(1) = 2.0;
      gp(2) = 0.0; 
	
    elseif points == 2	% 2-point Gauss integration in 1D

      gp     = zeros(2,2);

      gp(1,1) =  1.0;
      gp(2,1) = -1/sqrt(3);

      gp(1,2) =  1.0;
      gp(2,2) =  1/sqrt(3);


    else

      disp('Integration scheme not programmed')

    end  

  elseif dim == 2

    if points == 1	% 1-point Gauss integration in 2D for triangle

      gp     = zeros(3,1);

      gp(1,1) = 0.5;                          % weight
      gp(2,1) = 1.0/3.0; gp(3,1) = 1.0/3.0;   % position

    elseif points == 3	% 3-point Gauss integration in 2D for triangle

      gp     = zeros(3,3);

      gp(1,1) = 1.0/6.0;
      gp(2,1) = 0.0; gp(3,1) = 0.5;

      gp(1,2) = 1.0/6.0;
      gp(2,2) = 0.5; gp(3,2) = 0.0;

      gp(1,3) = 1.0/6.0;
      gp(2,3) = 0.5; gp(3,3) = 0.5;

    elseif points == 4	% 2x2-point Gauss integration in 2D for quad

      gp     = zeros(3,4);

      gp(1,:) = 1.0;

      gp(2,1) = -1.0/sqrt(3.0); gp(3,1) = -1.0/sqrt(3.0); 
      gp(2,2) =  1.0/sqrt(3.0); gp(3,2) = -1.0/sqrt(3.0); 
      gp(2,3) =  1.0/sqrt(3.0); gp(3,3) =  1.0/sqrt(3.0); 
      gp(2,4) = -1.0/sqrt(3.0); gp(3,4) =  1.0/sqrt(3.0); 

    elseif points == 9 % 3x3-point Gauss integration in 2D for quad

      gp     = zeros(3,9);
      a = sqrt(3.0/5.0);

      gp(1,1) = 25.0/81.0;
      gp(2,1) = -a; gp(3,1) = -a;

      gp(1,2) = 25.0/81.0;
      gp(2,2) = a; gp(3,2) = -a;

      gp(1,3) = 25.0/81.0;
      gp(2,3) = a; gp(3,3) = a;

      gp(1,4) = 25.0/81.0;
      gp(2,4) = -a; gp(3,4) = a;

      gp(1,5) = 40.0/81.0;
      gp(2,5) = 0.0; gp(3,5) = -a;

      gp(1,6) = 40.0/81.0;
      gp(2,6) = a; gp(3,6) = 0.0;

      gp(1,7) = 40.0/81.0;
      gp(2,7) = 0.0; gp(3,7) = a;

      gp(1,8) = 40.0/81.0;
      gp(2,8) = -a; gp(3,8) = 0.0;

      gp(1,9) = 64.0/81.0;
      gp(2,9) = 0.0; gp(3,9) = 0.0;

    else

      disp('Integration scheme not programmed')

    end

  elseif dim == 3

    disp('Integration scheme not programmed in 3D')

  end	