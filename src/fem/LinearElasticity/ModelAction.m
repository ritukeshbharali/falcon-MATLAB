%-------------------------------------------------------------------------%
% ModelAction is a function that carries out an action requested by the 
% solver. Here is a list of actions for the Linear Elasticity FE Model:
%    - INIT
%    - GET_MATRIX0
%    - GET_INT_FORCE
%    - GET_MATRIX2
%    - COMMIT
%    - REVERT
%    - REFINE_N_TRANSFER
%-------------------------------------------------------------------------%

function globdat = ModelAction(props,globdat,action)

switch action

    case 'INIT'

        globdat = init(props,globdat);

    case 'GET_MATRIX0'
        
        [globdat.K, globdat.fint] = ...
            getMatrix0(props,globdat);    

    case 'GET_INT_FORCE'
        
        globdat.fint = ...
            getIntForce(props,globdat);     
        
    case 'GET_MATRIX20'
        
        globdat.M = ...
            getMatrix20(props,globdat);
        
    case 'GET_MATRIX21'
        
        globdat.M = ...
            getMatrix21(props,globdat);    
        
    case 'COMMIT'    
        
        globdat.state00 = globdat.state0;
        globdat.state0  = globdat.state;

    case 'REVERT'    
        
        globdat.state   = globdat.state0;
        globdat.state0  = globdat.state00;    

    case 'REFINE_N_TRANSFER'

        % Refine and solution transfer
        globdat = refineNTransferSolution(props,globdat);

        % Recompute constraint matrix
        globdat = getConstraints(props,globdat);        

end

end


%=========================================================================%
% INIT  
%=========================================================================%

function globdat = init(props,globdat)

% Idea: Precompute model specific things required repeatedly.

end


%=========================================================================%
% GET_MATRIX0
%=========================================================================%

function [K, fint] = getMatrix0(props,globdat)

% Some useful variables
dim      = globdat.dim;
totdofs  = globdat.mesh.info.nnodes * dim;
nelems   = globdat.mesh.info.nelems.(props.feModel.elType);
elNodes  = length(globdat.mesh.elems.(props.feModel.elType)(1,2:end));

% Compute Gauss point locations and weights
gp = integrationScheme(props.feModel.Ipoints, globdat.dim);

% Sparse entries count 
idxCountK = nelems * elNodes^2 * dim^2;
idxCountf = nelems * elNodes * dim;

% Allocate coordinates of matrix entries and internal force vector
II = zeros(1,idxCountK);
JJ = zeros(1,idxCountK);
S  = zeros(1,idxCountK);

IIf = zeros(1,idxCountf);
Sf  = zeros(1,idxCountf);

% Loop over elements
for iel = 1:nelems

    % Element connectivity
    elConnect  = globdat.mesh.elems.(props.feModel.elType)(iel,2:end);
    elCoords   = globdat.mesh.nodes(elConnect,2:end);
    elUDof     = getDofMap(elConnect,globdat.dim);
    elNodes    = length(elConnect);

    % Get element displacement
    elDisp     = globdat.state(elUDof);

    % Set local matrix to zero
    elemMat    = zeros(length(elConnect)*dim);
    elfint     = zeros(length(elConnect)*dim,1);

    % Loop over Gauss points
    for gaussPoint = gp

        % Compute shape functions
        [N, dN, j] = shapeFunction(elCoords',gaussPoint,globdat.elType);

        % Compute strain and stress
        [~,Bd] = getDimShapeGrads(N,dN,dim);
        strain = Bd * elDisp;
        stress = props.mat.D * strain;

        % Compute element matrix and internal force
        dX = gaussPoint(1) * j;

        elemMat  = elemMat + Bd' * props.mat.D * Bd * dX;
        elfint   = elfint  + Bd' * stress * dX;  

    end

    elUDof    = reshape(elUDof,1,dim * elNodes);

    idxCountK  = (iel-1)*(dim)^2 * (length(elConnect))^2;
    idxCountf  = (iel-1)* dim * length(elConnect);

    % Get indexes and create row, col, vals for sparse assembly
    II(idxCountK+1:idxCountK+dim^2*elNodes^2) = repmat(elUDof,1,dim*elNodes);
    JJ(idxCountK+1:idxCountK+dim^2*elNodes^2) = reshape(repmat(elUDof,dim*elNodes,1),1,dim^2*elNodes^2);
    S(idxCountK+1:idxCountK+dim^2*elNodes^2)  = reshape(elemMat,1,dim^2*elNodes^2);

    IIf(idxCountf+1:idxCountf+dim*elNodes) = elUDof;
    Sf(idxCountf+1:idxCountf+dim*elNodes)  = elfint;

end

% Assemble (sparse) stiffness matrix
K    = sparse(II,JJ,S,totdofs,totdofs);
fint = sparse(IIf,ones(length(IIf),1),Sf,totdofs,1);

end



%=========================================================================%
%  GET_INT_FORCE
%=========================================================================%

function [fint] = getIntForce(props,globdat)

fint        = [];

end



%=========================================================================%
% GET_MATRIX20
%=========================================================================%

function M = getMatrix20(props,globdat)

M       = [];

end


%=========================================================================%
% GET_MATRIX21
%=========================================================================%

function M = getMatrix21(props,globdat)

M       = [];

end

