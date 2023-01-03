% Copyright (C) 2007 Garth N. Wells
%
% Write VTK post-processing files


function VTKPostProcess(dir, props, globdat)

dim      = globdat.dim;
numNodes = globdat.mesh.info.nnodes;
numCells = globdat.mesh.info.nelems.(props.feModel.elType);
dof      = globdat.dim; 
etype    = 'Tri3';
x        = transpose(globdat.mesh.nodes(:,2:3));
connect  = transpose(globdat.mesh.elems.(props.feModel.elType)(:,2:4));
step     = globdat.ts.step;

% Output file extension
if( step < 10 )
  ext = '0000';
elseif( step >= 10 && step < 100 ) 
  ext = '000';
elseif( step >= 100 && step < 1000 ) 
  ext = '00';
else
  error('To many steps for VTK output')
end

% Output files
outfilePVD  = [dir '/results.pvd'];
outfileVTU  = [dir '/results' ext int2str(step) '.vtu'];
results_vtu = fopen(outfileVTU, 'wt');

fileVTU = ['results' ext int2str(step) '.vtu'];

if( step == 0 )
  results_pvd = fopen(outfilePVD, 'wt');
  fprintf(results_pvd, '<?xml version="1.0"?> \n'); 
  fprintf(results_pvd, '<VTKFile type="Collection" version="0.1" > \n'); 
  fprintf(results_pvd, '<Collection> \n'); 
  fprintf(results_pvd, '<DataSet timestep="%g" part="0" file="%s"/> \n', step, fileVTU); 
  fprintf(results_pvd, '</Collection> \n'); 
  fprintf(results_pvd, '</VTKFile> \n'); 
else

  results_pvd = fopen(outfilePVD, 'r+t');
  found = 0;
  while(~found)
    old_position = ftell(results_pvd);
    line = fgetl(results_pvd);
    found = length(findstr(line, '</Collection>'));
    new_position = ftell(results_pvd);
  end
  fseek(results_pvd, old_position, 'bof'); 

  fprintf(results_pvd, '<DataSet timestep="%g" part="0" file="%s"/> \n', step, fileVTU); 
  fprintf(results_pvd, '</Collection> \n'); 
  fprintf(results_pvd, '</VTKFile> \n'); 

end


if(strcmp(etype, 'Bar1') || strcmp(etype, 'Bar2')) 
  numVertexesPerCell = 2;
  VTKCellCode = 3;
elseif(strcmp(etype, 'Quad4') || strcmp(etype, 'Quad8') || strcmp(etype, 'Quad9')) 
  numVertexesPerCell = 4;
  VTKCellCode = 9;
elseif(strcmp(etype, 'Tri3') || strcmp(etype, 'Tri6')) 
  numVertexesPerCell = 3;
  VTKCellCode = 5;
elseif(strcmp(etype, 'Tet4'))
  numVertexesPerCell = 4;
  VTKCellCode = 10;
else
  error('Element type not known (VTKPostProcess)')
end

% dof_per_vertex = 0;
if (dof == 1)
  dof_per_vertex = 1;
elseif (dof == 2 )
  dof_per_vertex = 2;
elseif (dof == 3 )
  dof_per_vertex = 3;
else
  dof_per_vertex = 3;
  disp('Only outputting first three degrees of freedom')
end

% Write headers
fprintf(results_vtu, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n'); 
fprintf(results_vtu, '<UnstructuredGrid> \n'); 
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfCells=" %g"> \n', numNodes, numCells); 

% Write point data
fprintf(results_vtu, '<Points> \n'); 
if( dof_per_vertex == 1)
  fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="1"  format="ascii" > \n'); 
else
  fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n'); 
end

  for i=1:numNodes
    if( dim == 3)
      fprintf(results_vtu, '%f ',  x(1:3,i));
    elseif(dim == 2)
      fprintf(results_vtu, '%f ',  x(1:2,i));
      fprintf(results_vtu, '0.0 ');
    elseif(dim == 1)
      fprintf(results_vtu, '%f ',  x(1:1,i));
      fprintf(results_vtu, '0.0 ');
      fprintf(results_vtu, '0.0 ');
    end
    fprintf(results_vtu, '\n');
  end						 

fprintf(results_vtu, '</DataArray> \n'); 
fprintf(results_vtu, '</Points> \n'); 

% Print cells 
fprintf(results_vtu, '<Cells> \n'); 

% Print cell connectivity
fprintf(results_vtu, '<DataArray  type="Int32"  Name="connectivity"  format="ascii"> \n'); 

  for i=1:numCells
    fprintf(results_vtu, '%g ',  connect(1:numVertexesPerCell, i)-1 );
  	fprintf(results_vtu, '\n');
  end						 

fprintf(results_vtu, '</DataArray> \n'); 

% Print cell offsets
fprintf(results_vtu, '<DataArray  type="Int32"  Name="offsets"  format="ascii"> \n'); 

  offset = 0;
  for i=1:numCells
    offset = offset + numVertexesPerCell;
    fprintf(results_vtu, '%g ', offset);
  	fprintf(results_vtu, '\n');
  end						 

fprintf(results_vtu, '</DataArray> \n'); 

% Print cell types
fprintf(results_vtu, '<DataArray  type="UInt8"  Name="types"  format="ascii"> \n'); 

  for i=1:numCells
    fprintf(results_vtu, '%g ', VTKCellCode);
  	fprintf(results_vtu, '\n');
  end

fprintf(results_vtu, '</DataArray> \n'); 
fprintf(results_vtu, '</Cells> \n');

% Split solution based on femodel and dimension
switch props.feModel.type

    case 'PhaseFieldFracture'

        u   = globdat.state(1:dim*numNodes);
        phi = globdat.state(dim*numNodes+1:(dim+1)*numNodes);

        % Print displacement
        fprintf(results_vtu, '<PointData  Vectors="U"> \n'); 
        fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n'); 

        for i=1:numNodes	
            if( dof_per_vertex ~= 2)
                for j=1:dof_per_vertex
                    fprintf(results_vtu, '%g   ', u( dof*(i-1) + j) );
                end
            else
                for j=1:dof_per_vertex
                    fprintf(results_vtu, '%g   ', u( dof*(i-1) + j) );
                end
            	fprintf(results_vtu, '  0.0');
            end
        fprintf(results_vtu, '\n');
        end
    
        fprintf(results_vtu, '</DataArray> \n'); 

        % Print phase-field 
        fprintf(results_vtu, '<DataArray  type="Float64"  Name="Phi" NumberOfComponents="1" format="ascii"> \n');

        for i=1:numNodes	
            fprintf(results_vtu, '%g   ', phi(i) );
            fprintf(results_vtu, '\n');
        end

        fprintf(results_vtu, '</DataArray> \n'); 
        fprintf(results_vtu, '</PointData> \n'); 

    otherwise
        error('Not yet implemented!')
end 

fprintf(results_vtu, '</Piece> \n'); 
fprintf(results_vtu, '</UnstructuredGrid> \n'); 
fprintf(results_vtu, '</VTKFile> \n'); 

fclose(results_vtu);
fclose(results_pvd);