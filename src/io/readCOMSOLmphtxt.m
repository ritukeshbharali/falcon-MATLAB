% ------------------------------------------------------------------------%
% FUNCTION PURPOSE: Extract mesh information from a COMSOL mesh file      % 
%                                                                         %
%    INPUT:  file_name - must have extension '.mphtxt'                    %
%                                                                         %
%    OUTPUT: nodes - contains node numbers and coordinates                %
%            elems - struct containing different element types with their %
%                    connectivity data                                    %
%            idx   - struct containing additional information tags        %
%                                                                         %
% Author: Ritukesh Bharali (ritukesh.bharali@chalmers.se)                 %
% Date:   19.11.2021                                                      %
% ----------------------------------------------------------------------- %

function [Mesh] = readCOMSOLmphtxt(file_name)

Mesh = struct;
info = struct;

file = file_name;
file_input   = fopen(file,'rt');

disp(' - Reading input file')

% ------------- Test input file ----------------------

% Read 2 first lines and check if we have correct format
line        = fgetl(file_input);
mesh_format = split(line);
comp        = 'COMSOL';

if strcmp (mesh_format{4},comp) == 1
    disp('  -- COMSOL mesh: Type mphtxt!')
else
    error('  -- Wrong file type, must be COMSOL .mphtxt file')
end


% ------------- Extract nodes ----------------------

disp('  -- Extracting nodal data ...')

while ~feof(file_input) 
    
    line  = fgetl(file_input);
    
    if isempty(line) == 1
        continue;
    end
    
    buf   = split(line);
     
     if length(buf) == 3 && strcmp (buf{3},'sdim') == 1
        disp(['       - Dimension: ',buf{1}])
        ndim      = str2double(buf{1});
        info.ndim = ndim;
     elseif length(buf) == 6 && strcmp (buf{6},'vertices') == 1
        disp(['       - No. of nodes: ',buf{1}])
        nnodes      = str2double(buf{1});
        info.nnodes = nnodes;
     elseif length(buf) == 4 && strcmp (buf{4},'coordinates') == 1
        break;
     end
end

% Gather all nodes and coordinates
nodes = gather_entity(file_input,nnodes,ndim+1,'coordinates');



% ------------- Extract elements ----------------------

disp('  -- Extracting element data ...')

while ~feof(file_input) 
    
    line  = fgetl(file_input);
    
    if isempty(line) == 1
        continue;
    end
    
    buf   = split(line);
    
    if length(buf) == 5 && strcmp (buf{3},'Object') == 1
        break;
    end
    
    if length(buf) == 3 && strcmp (buf{2},'Type') == 1
        
        % Get element type
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        buf   = split(line);
        name  = buf{2};
        disp(['      Found element type: ',name])
        
        % Get element vertices
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        buf   = split(line);
        vert  = str2double(buf{1}); % Vertices
        
        % Get number of elements
        line  = fgetl(file_input);
        buf   = split(line);
        nelem = str2double(buf{1}); 
        
        % Extract element data
        line  = fgetl(file_input);
        el    = gather_entity(file_input,nelem,vert+1,'nodes');
        varname = matlab.lang.makeValidName(name);
        elems.(varname) = el;
        
        % Extract associated labels
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        idxel = gather_index(file_input,nelem);
        idx.(varname) = idxel;
        
        info.nelems.(varname) = nelem;
        
        disp(['       - No. of elements: ',num2str(nelem)])
        
    end
    
end


   
% ------------- Extract additonal labels ----------------------

disp('  -- Extracting additional labels ...')

while ~feof(file_input)
    
    line  = fgetl(file_input);
    
    if isempty(line) == 1
        continue;
    end
    
    buf = split(line);
    
    if length(buf) == 4 && strcmp (buf{2},'Selection') == 1 
        
        % Get what the selection is
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        buf   = split(line);
        name  = buf{2};
        disp(['      Found label: ',name])
        
        % Extract size
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        line  = fgetl(file_input);
        buf   = split(line);
        n     = str2double(buf{1});
        line  = fgetl(file_input);
        
        % Extract data
        data = zeros(n,1);
        for i = 1:n
            line       = fgetl(file_input);
            buf        = split(line);
            data(i)    = str2double(buf{1});
        end
        
        varname = matlab.lang.makeValidName(name);
        idx.(varname) = data;
        
    end
    
    
end
  
disp(' - Completed mesh extraction')

Mesh.nodes = nodes;
Mesh.elems = elems;
Mesh.idx   = idx;
Mesh.info  = info;

end









% --- Function to gather entities (nodes, elements) ---- %

function entity = gather_entity(file_input,nrows,ncols,type)

entity = zeros(nrows,ncols);

for i = 1:nrows

line  = fgetl(file_input);
buf   = split(line);
    
% Fill first col
entity(i,1) = i;

    % Fill remaining cols
    for j = 1:ncols-1
        
        switch type
           
            case 'none'
                entity(i,j+1) = str2double(buf{j});
                
            case 'coordinates'
                entity(i,j+1) = str2double(buf{j});
                
            case 'nodes'
                
                % + 1 since COMSOL node numbering starts from one
                entity(i,j+1) = str2double(buf{j}) + 1;
                
            otherwise
                error('Pick coordinates or elements!')
        end
        
    end
end

end


% --- Function to gather index (index lists) ---- %

function entity = gather_index(file_input,nrows)

entity = zeros(nrows,1);

for i = 1:nrows

line  = fgetl(file_input);
buf   = split(line);
    
% Fill first col
entity(i,1) = str2double(buf{1});

end

end