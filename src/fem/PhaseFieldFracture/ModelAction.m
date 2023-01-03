%-------------------------------------------------------------------------%
% ModelAction is a function that carries out an action requested by the 
% solver. Here is a list of actions for the Phase-Field Fracture FE Model:
%    - INIT
%    - GET_MATRIX0
%    - GET_MATRIX0_EXT
%    - GET_INT_FORCE
%    - GET_MATRIX0_1
%    - GET_MATRIX0_2
%    - GET_DISSIPATION
%    - GET_AT2_DISSIPATION
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

    case 'GET_MATRIX0_EXT'

        [globdat.K, globdat.fint] = ...
            getMatrix0_Ext(props,globdat);     

    case 'GET_INT_FORCE'
        
        globdat.fint = ...
            getIntForce(props,globdat);     
    
    case 'GET_MATRIX0_1'
        
        [globdat.K1, globdat.fint1] = ...
            getMatrix0_1(props,globdat);
        
    case 'GET_MATRIX0_2'
        
        [globdat.K2, globdat.fint2] = ...
            getMatrix0_2(props,globdat);

    case 'GET_DISSIPATION'
        
        [globdat.h, globdat.g] = ...
            getDissipation(props,globdat);

    case 'GET_AT2_DISSIPATION'
        
        [globdat.h, globdat.g] = ...
            getAT2Dissipation(props,globdat);    
        
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
totdofs  = globdat.mesh.info.nnodes * (dim+1);
udofs    = globdat.mesh.info.nnodes * dim;
nelems   = globdat.mesh.info.nelems.(props.feModel.elType);
elNodes  = length(globdat.mesh.elems.(props.feModel.elType)(1,2:end));

% Compute Gauss point locations and weights
gp = integrationScheme(props.feModel.Ipoints, globdat.dim);

% Sparse entries count 
idxCountK = nelems * elNodes^2 * (dim+1)^2;
idxCountf = nelems * elNodes * (dim+1);

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
    elPhiDof   = elConnect + udofs;
    elNodes    = length(elConnect);

    % Get element phase-field and displacement
    elPhi      = globdat.state(elPhiDof);
    elDisp     = globdat.state(elUDof);

    % Set local matrix to zero
    elemMat11   = zeros(length(elConnect)*dim);
    elemMat12   = zeros(length(elConnect)*dim,length(elConnect));
    elemMat21   = zeros(length(elConnect),length(elConnect)*dim);
    elemMat22   = zeros(length(elConnect));

    elfint1     = zeros(length(elConnect)*dim,1);
    elfint2     = zeros(length(elConnect),1);

    % Loop over Gauss points
    for gaussPoint = gp

        % Compute shape functions
        [N, dN, j] = shapeFunction(elCoords',gaussPoint,globdat.elType);

        % Compute phase-field at Gauss point and compute degradation
        gp_phi     = N * elPhi;
        grad_phi   = dN * elPhi;

        % Compute degradation and its derivatives
        deg_phi    = (1-gp_phi)^2;
        dg_phi   = -2.0 * (1.0 - gp_phi);
        ddg_phi  = 2.0;

        % Compute strain and stress
        [~,Bd] = getDimShapeGrads(N,dN,dim);
        strain = Bd * elDisp;
        stress = props.mat.D * strain;

        switch props.mat.Esplit
            
            case 'NoSplit'
                
                Psi = 0.5 * dot(stress,strain);
                
            case 'Spectral2D'  % Miehe et. al. (2010)
                strain_tensor = [strain(1) strain(3)/2; strain(3)/2 strain(2)];  
                tracep_strain = max(0,strain(1)+strain(2));
                eigp_strain   = max(0,eig(strain_tensor));
                Psi           = 0.5 * props.mat.lambda * tracep_strain^2 + ...
                                props.mat.mu*(eigp_strain(1)^2 + eigp_strain(2)^2); 

            case 'Amor'        % Amor et. al. (2009)
                tracep_strain = max(0,strain(1)+strain(2));
                vol_strain    = 0.5*(strain(1)+strain(2));
                dev_strain    = [strain(1)-vol_strain;
                                 strain(2)-vol_strain;
                                 strain(3)/2];
                Psi           = 0.5 * props.mat.K * tracep_strain^2 + ...
                                props.mat.mu * dot(dev_strain,dev_strain);
                
            otherwise
                error('Not implemented!')
        end

        % Compute element matrix and internal force
        dX = gaussPoint(1) * j;

        elemMat11  = elemMat11 + deg_phi * Bd' * props.mat.D * Bd * dX;
        elemMat12  = elemMat12 + dg_phi  * Bd' * stress * N * dX;
        elemMat21  = elemMat21 + dg_phi  * N' *  stress' * Bd * dX;
        elemMat22  = props.mat.gc * props.mat.l0 ...
                    * (dN'*dN) * dX;
        elemMat22  = elemMat22 + (props.mat.gc/props.mat.l0 + ...
                      ddg_phi * Psi ) * (N'*N) * dX;

        elfint1    = elfint1 + deg_phi * Bd' * stress * dX;  
        elfint2    = elfint2 + props.mat.gc*props.mat.l0 * ...
                      (dN' * grad_phi ) * dX;
        elfint2    = elfint2 + (props.mat.gc/props.mat.l0*gp_phi + ...
                      dg_phi * Psi ) * N' * dX;

    end

    idxCountK  = (iel-1)*(dim+1)^2 * (length(elConnect))^2;
    idxCountf  = (iel-1)* (dim+1) * length(elConnect);

    sctrx     = elConnect;
    dsctrx    = [reshape([2*sctrx-1;2*sctrx],1,(dim)*elNodes), udofs + sctrx];

    % Get indexes and create row, col, vals for sparse assembly
    II(idxCountK+1:idxCountK+(dim+1)^2*elNodes^2) = repmat(dsctrx,1,(dim+1)*elNodes);
    JJ(idxCountK+1:idxCountK+(dim+1)^2*elNodes^2) = reshape(repmat(dsctrx,(dim+1)*elNodes,1),1,(dim+1)^2*elNodes^2);
    S(idxCountK+1:idxCountK+(dim+1)^2*elNodes^2)  = reshape([elemMat11 elemMat12; elemMat21 elemMat22],1,(dim+1)^2*elNodes^2);

    IIf(idxCountf+1:idxCountf+(dim+1)*elNodes) = dsctrx;
    Sf(idxCountf+1:idxCountf+(dim+1)*elNodes)  = [elfint1;elfint2];

end

% Assemble (sparse) stiffness matrix
K    = sparse(II,JJ,S,totdofs,totdofs);
fint = sparse(IIf,ones(length(IIf),1),Sf,totdofs,1);

end



%=========================================================================%
%  GET_MATRIX0_EXT
%=========================================================================%

function [K, fint] = getMatrix0_Ext(props,globdat)

% Some useful variables
dim      = globdat.dim;
totdofs  = globdat.mesh.info.nnodes * (dim+1);
udofs    = globdat.mesh.info.nnodes * dim;
nelems   = globdat.mesh.info.nelems.(props.feModel.elType);
elNodes  = length(globdat.mesh.elems.(props.feModel.elType)(1,2:end));

% Compute Gauss point locations and weights
gp = integrationScheme(props.feModel.Ipoints, globdat.dim);

% Sparse entries count 
idxCountK = nelems * elNodes^2 * (dim+1)^2;
idxCountf = nelems * elNodes * (dim+1);

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
    elPhiDof   = elConnect + udofs;
    elNodes    = length(elConnect);

    % Get element phase-field and displacement
    elPhi      = globdat.state(elPhiDof);
    elPhi0     = globdat.state0(elPhiDof);
    elPhiEx    = globdat.stateEx(elPhiDof);
    elDisp     = globdat.state(elUDof);

    % Set local matrix to zero
    elemMat11   = zeros(length(elConnect)*dim);
    elemMat12   = zeros(length(elConnect)*dim,length(elConnect));
    elemMat21   = zeros(length(elConnect),length(elConnect)*dim);
    elemMat22   = zeros(length(elConnect));

    elfint1     = zeros(length(elConnect)*dim,1);
    elfint2     = zeros(length(elConnect),1);

    % Loop over Gauss points
    for gaussPoint = gp

        % Compute shape functions
        [N, dN, j] = shapeFunction(elCoords',gaussPoint,globdat.elType);

        % Compute phase-field at Gauss point and compute degradation
        gp_phi     = N * elPhi;
        gp_phi0    = N * elPhi0;
        gp_phiEx   = N * elPhiEx;
        grad_phi   = dN * elPhi;

        % Compute degradation and its derivatives
        deg_phiEx  = (1-gp_phiEx)^2;
        dg_phi     = -2.0 * (1.0 - gp_phi);
        ddg_phi    = 2.0;

        % Compute strain and stress
        [~,Bd] = getDimShapeGrads(N,dN,dim);
        strain = Bd * elDisp;
        stress = props.mat.D * strain;

        switch props.mat.Esplit
            
            case 'NoSplit'
                
                Psi = 0.5 * dot(stress,strain);
                
            case 'Spectral2D'  % Miehe et. al. (2010)
                strain_tensor = [strain(1) strain(3)/2; strain(3)/2 strain(2)];  
                tracep_strain = max(0,strain(1)+strain(2));
                eigp_strain   = max(0,eig(strain_tensor));
                Psi           = 0.5 * props.mat.lambda * tracep_strain^2 + ...
                                props.mat.mu*(eigp_strain(1)^2 + eigp_strain(2)^2); 

            case 'Amor'        % Amor et. al. (2009)
                tracep_strain = max(0,strain(1)+strain(2));
                vol_strain    = 0.5*(strain(1)+strain(2));
                dev_strain    = [strain(1)-vol_strain;
                                 strain(2)-vol_strain;
                                 strain(3)/2];
                Psi           = 0.5 * props.mat.K * tracep_strain^2 + ...
                                props.mat.mu * dot(dev_strain,dev_strain);
                
            otherwise
                error('Not implemented!')
        end

        % Compute element matrix and internal force
        dX = gaussPoint(1) * j;

        elemMat11  = elemMat11 + deg_phiEx * Bd' * props.mat.D * Bd * dX;
        elemMat21  = elemMat21 + dg_phi  * N' *  stress' * Bd * dX;
        elemMat22  = props.mat.gc * props.mat.l0 ...
                    * (dN'*dN) * dX;
        elemMat22  = elemMat22 + (props.mat.gc/props.mat.l0 + ...
                      ddg_phi * Psi ) * (N'*N) * dX;

        elfint1    = elfint1 + deg_phiEx * Bd' * stress * dX;  
        elfint2    = elfint2 + props.mat.gc*props.mat.l0 * ...
                      (dN' * grad_phi ) * dX;
        elfint2    = elfint2 + (props.mat.gc/props.mat.l0*gp_phi + ...
                      dg_phi * Psi ) * N' * dX;

    end

    idxCountK  = (iel-1)*(dim+1)^2 * (length(elConnect))^2;
    idxCountf  = (iel-1)* (dim+1) * length(elConnect);

    sctrx     = elConnect;
    dsctrx    = [reshape([2*sctrx-1;2*sctrx],1,(dim)*elNodes), udofs + sctrx];

    % Get indexes and create row, col, vals for sparse assembly
    II(idxCountK+1:idxCountK+(dim+1)^2*elNodes^2) = repmat(dsctrx,1,(dim+1)*elNodes);
    JJ(idxCountK+1:idxCountK+(dim+1)^2*elNodes^2) = reshape(repmat(dsctrx,(dim+1)*elNodes,1),1,(dim+1)^2*elNodes^2);
    S(idxCountK+1:idxCountK+(dim+1)^2*elNodes^2)  = reshape([elemMat11 elemMat12; elemMat21 elemMat22],1,(dim+1)^2*elNodes^2);

    IIf(idxCountf+1:idxCountf+(dim+1)*elNodes) = dsctrx;
    Sf(idxCountf+1:idxCountf+(dim+1)*elNodes)  = [elfint1;elfint2];

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
% GET_MATRIX0_1
%=========================================================================%

function [K, fint] = getMatrix0_1(props,globdat)

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
    elPhiDof   = elConnect + globdat.cons.nDof1;
    elNodes    = length(elConnect);

    % Get element phase-field and displacement
    elPhi      = globdat.state(elPhiDof);
    elDisp     = globdat.state(elUDof);

    % Set local matrix to zero
    elemMat    = zeros(length(elConnect)*dim);
    elfint = zeros(length(elConnect)*dim,1);

    % Loop over Gauss points
    for gaussPoint = gp

        % Compute shape functions
        [N, dN, j] = shapeFunction(elCoords',gaussPoint,globdat.elType);

        % Compute phase-field at Gauss point and compute degradation
        gp_phi     = N * elPhi;
        deg_phi    = (1-gp_phi)^2;

        % Compute strain and stress
        [~,Bd] = getDimShapeGrads(N,dN,dim);
        strain = Bd * elDisp;
        stress = props.mat.D * strain;

        % Compute element matrix and internal force
        dX = gaussPoint(1) * j;
        elemMat    = elemMat    + deg_phi * Bd' * props.mat.D * Bd * dX;
        elfint = elfint + deg_phi * Bd' * stress * dX;  

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
% GET_MATRIX0_2
%=========================================================================%

function [K, fint] = getMatrix0_2(props,globdat)

% Some useful variables
dim      = globdat.dim;
totdofs  = globdat.mesh.info.nnodes;
nelems   = globdat.mesh.info.nelems.(props.feModel.elType);
elnodes  = length(globdat.mesh.elems.(props.feModel.elType)(1,2:end));

% Compute Gauss point locations and weights
gp = integrationScheme(props.feModel.Ipoints, dim);

% Sparse entries count 
idxCountK = nelems * elnodes^2;
idxCountf = nelems * elnodes;

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
    elUDof     = getDofMap(elConnect,dim);
    elPhiDof   = elConnect + globdat.cons.nDof1;

    % Get element phase-field and displacement
    elPhi      = globdat.state(elPhiDof);
    elDisp     = globdat.state(elUDof);

    % Set local matrix to zero
    elemMat    = zeros(length(elConnect));
    elfint = zeros(length(elConnect),1);

    % Loop over Gauss points
    for gaussPoint = gp
        % gauss_point = gp;

        % Compute shape functions
        [N, dN, j] = shapeFunction(elCoords',gaussPoint,globdat.elType);

        % Compute phase-field at Gauss point and compute degradation
        gp_phi     = N * elPhi;
        grad_phi   = dN * elPhi;

        % Compute degradation and its derivatives
        dg_phi   = -2.0 * (1.0 - gp_phi);
        ddg_phi  = 2.0;

        % Compute strain and stress
        [~,Bd] = getDimShapeGrads(N,dN,dim);
        strain = Bd * elDisp;
        stress = props.mat.D * strain;

        switch props.mat.Esplit
            
            case 'NoSplit'
                
                Psi = 0.5 * dot(stress,strain);
                
            case 'Spectral2D'  % Miehe et. al. (2010)
                strain_tensor = [strain(1) strain(3)/2; strain(3)/2 strain(2)];  
                tracep_strain = max(0,strain(1)+strain(2));
                eigp_strain   = max(0,eig(strain_tensor));
                Psi           = 0.5 * props.mat.lambda * tracep_strain^2 + ...
                                props.mat.mu*(eigp_strain(1)^2 + eigp_strain(2)^2); 

            case 'Amor'        % Amor et. al. (2009)
                tracep_strain = max(0,strain(1)+strain(2));
                vol_strain    = 0.5*(strain(1)+strain(2));
                dev_strain    = [strain(1)-vol_strain;
                                 strain(2)-vol_strain;
                                 strain(3)/2];
                Psi           = 0.5 * props.mat.K * tracep_strain^2 + ...
                                props.mat.mu * dot(dev_strain,dev_strain);
                
            otherwise
                error('Not implemented!')
        end

        % Compute element matrix and internal force
        dX = gaussPoint(1) * j;
        elemMat    = elemMat + props.mat.gc * props.mat.l0 ...
                    * (dN'*dN) * dX;
        elemMat    = elemMat + (props.mat.gc/props.mat.l0 + ...
                      ddg_phi * Psi ) * (N'*N) * dX;
        elfint = elfint + props.mat.gc*props.mat.l0 * ...
                      (dN' * grad_phi ) * dX;
        elfint = elfint + (props.mat.gc/props.mat.l0*gp_phi + ...
                      dg_phi * Psi ) * N' * dX;
    end

    idxCountK  = (iel-1)* (length(elConnect))^2;
    idxCountf  = (iel-1)* length(elConnect);

    % Get indexes and create row, col, vals for sparse assembly
    II(idxCountK+1:idxCountK+elnodes^2) = repmat(elConnect,1,elnodes);
    JJ(idxCountK+1:idxCountK+elnodes^2) = reshape(repmat(elConnect,elnodes,1),1,elnodes^2);
    S(idxCountK+1:idxCountK+elnodes^2)  = reshape(elemMat,1,elnodes^2);

    IIf(idxCountf+1:idxCountf+elnodes) = elConnect;
    Sf(idxCountf+1:idxCountf+elnodes)  = elfint;

end

% Assemble (sparse) stiffness matrix
K    = sparse(II,JJ,S,totdofs,totdofs);
fint = sparse(IIf,ones(length(IIf),1),Sf,totdofs,1);

end


%=========================================================================%
% GET_DISSIPATION
%=========================================================================%

function [h, g] = getDissipation(props,globdat)

% Initialise h
h = 0;
g = 0;

end




%=========================================================================%
% GET_AT2_DISSIPATION
%=========================================================================%

function [h, g] = getAT2Dissipation(props,globdat)

% Initialise h
h = 0;
g = 0;

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

