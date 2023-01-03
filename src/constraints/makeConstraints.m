function [cons,gdofs] = makeConstraints(props,globdat)

% Check if model has constraints
ncons = length(props.cons);
assert(ncons>0,...
    'Without constraints, the system of equations would be singular!')

% Allocate some data
gdofs = [];
gvals = [];
ldofs = [];

% Get FE model type
switch props.feModel.type

    case 'Elasticity'

        for i = 1:ncons
            
            % Throw error if a particular edge label does not exist
            assert(isfield(globdat.mesh.idx,props.cons(i).edge),...
            'The chosen edge is not found in the imported mesh!')

            % Extract constrained nodes
            cons_el    = find(ismember(globdat.mesh.idx.edg, ...
                         globdat.mesh.idx.(props.cons(i).edge{1,1})) == 1);
            cons_nodes = globdat.mesh.elems.edg(cons_el,2:end);
            cons_nodes = unique(cons_nodes);

           if isrow(cons_nodes)
               cons_nodes = cons_nodes';
           end

           nConsDofs = length(props.cons(i).dofs);
           assert(nConsDofs>0,...
           'Without constraints, the system of equations would be singular!')

           for j = 1:nConsDofs

               dof_type = props.cons(i).dofs{1,j};
               vals     = props.cons(i).vals{1,j};

               switch globdat.dim

                   case 2

                       if strcmp(dof_type,'dx') == 1

                           gdofs = [gdofs; 2*cons_nodes-1];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 2*cons_nodes-1];
                           end

                       elseif strcmp(dof_type,'dy') == 1

                           gdofs = [gdofs; 2*cons_nodes];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 2*cons_nodes];
                           end

                       end

                   case 3

                       if strcmp(dof_type,'dx') == 1

                           gdofs = [globdat.cons.dofs; 3*cons_nodes-2];
                           gvals = [globdat.cons.vals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 3*cons_nodes-2];
                           end

                       elseif strcmp(dof_type,'dy') == 1

                           gdofs = [gdofs; 3*cons_nodes-1];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 3*cons_nodes-1];
                           end

                       elseif strcmp(dof_type,'dz') == 1

                           gdofs = [gdofs; 3*cons_nodes];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];    

                           if vals == 1
                               ldofs = [ldofs; 3*cons_nodes];
                           end
                       end

                   otherwise
                       error('Fix bug in code, error should be thrown earlier!')
               end
           end
        end

    case 'PhaseFieldFracture'

        for i = 1:ncons
            
            % Throw error if a particular edge label does not exist
            assert(isfield(globdat.mesh.idx,props.cons(i).edge),...
            'The chosen edge is not found in the imported mesh!')

            % Extract constrained nodes
            cons_el    = find(ismember(globdat.mesh.idx.edg, ...
                         globdat.mesh.idx.(props.cons(i).edge{1,1})) == 1);
            cons_nodes = globdat.mesh.elems.edg(cons_el,2:end);
            cons_nodes = unique(cons_nodes);

           if isrow(cons_nodes)
               cons_nodes = cons_nodes';
           end

           nConsDofs = length(props.cons(i).dofs);
           assert(nConsDofs>0,...
           'Without constraints, the system of equations would be singular!')

           for j = 1:nConsDofs

               dof_type = props.cons(i).dofs{1,j};
               vals     = props.cons(i).vals{1,j};

               switch globdat.dim

                   case 2

                       if strcmp(dof_type,'dx') == 1

                           gdofs = [gdofs; 2*cons_nodes-1];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 2*cons_nodes-1];
                           end

                       elseif strcmp(dof_type,'dy') == 1

                           gdofs = [gdofs; 2*cons_nodes];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 2*cons_nodes];
                           end

                       end

                   case 3

                       if strcmp(dof_type,'dx') == 1

                           gdofs = [globdat.cons.dofs; 3*cons_nodes-2];
                           gvals = [globdat.cons.vals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 3*cons_nodes-2];
                           end

                       elseif strcmp(dof_type,'dy') == 1

                           gdofs = [gdofs; 3*cons_nodes-1];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];

                           if vals == 1
                               ldofs = [ldofs; 3*cons_nodes-1];
                           end

                       elseif strcmp(dof_type,'dz') == 1

                           gdofs = [gdofs; 3*cons_nodes];
                           gvals = [gvals; vals*ones(length(cons_nodes),1)];    

                           if vals == 1
                               ldofs = [ldofs; 3*cons_nodes];
                           end
                       end

                   otherwise
                       error('Fix bug in code, error should be thrown earlier!')
               end
           end
        end

    otherwise
        error('BUG: Should have thrown a wrong FE model error earlier!')

end

% Store constrained dofs and values in cons in globdat
cons.tab   = [gdofs, gvals];

% Store loaded dofs for lodi plots
cons.ldofs = ldofs;

end