function globdat = initProblem(props)

% Create the globdat struct
globdat = struct;

% Import FE mesh
globdat.mesh = readCOMSOLmphtxt(props.meshFile);

% Color mesh
%[globdat.mesh.info.C,globdat.mesh.info.nColors] = GreedyMeshColoring(...
%    globdat.mesh.elems.(props.feModel.elType)(:,2:end),8);

% Get FE model type
switch props.feModel.type

    case 'Elasticity'
        globdat.ndofs = globdat.mesh.info.ndim;
        globdat.dim   = globdat.mesh.info.ndim;

        % Initialise vector for load-displacement plot
        globdat.lodi = zeros(props.ts.nSteps+1,2);

    case 'PhaseFieldFracture'
        globdat.ndofs = globdat.mesh.info.ndim + 1;
        globdat.dim   = globdat.mesh.info.ndim;

        % Initialise vector for load-displacement plot
        globdat.lodi = zeros(props.ts.nSteps+1,2);

    otherwise
        error('Wrong choice of FE Model!')

end

% Ask FE model to initialize
globdat = ModelAction(props,globdat,'INIT');

% Check if problem element type exists in the mesh
assert(isfield(globdat.mesh.elems,props.feModel.elType),...
    'The chosen element type is not found in the imported mesh!')

% Get element type
elNumNodes     = length(globdat.mesh.elems.(props.feModel.elType)(1,2:end));
globdat.elType = append(props.feModel.elType,num2str(elNumNodes));
globdat.elType = replace(globdat.elType,globdat.elType(1),upper(globdat.elType(1)));

% Allocate current, old step. old old step solution vectors
ntotdofs         = globdat.ndofs * globdat.mesh.info.nnodes;
globdat.state    = zeros(ntotdofs,1);
globdat.state0   = zeros(ntotdofs,1);
globdat.state00  = zeros(ntotdofs,1);

switch props.nlSolver.type

    case 'newtonEx'
        globdat.stateEx = zeros(ntotdofs,1);
        globdat.ts.dt0  = 0.0;
        globdat.ts.dt00 = 0.0;
end

% Check if model has constraints and make constraints in globdat
[globdat.cons,gdofs] = makeConstraints(props,globdat);

% Initialize some parameters for timestepping
globdat.ts.step = 0;
globdat.ts.t    = 0.0;
globdat.ts.dt   = props.ts.initStepSize;

% Set globdat status
globdat.active = true;
globdat.redo   = false;

% Create a constraint matrix (to be removed later)
if strcmp(props.nlSolver.type,'staggered')
    globdat.cons.nDof1 = globdat.dim * globdat.mesh.info.nnodes;
    globdat.cons.nDof2 = (globdat.ndofs-globdat.dim) * ...
                          globdat.mesh.info.nnodes;
    globdat.cons.CMat1 = speye(globdat.cons.nDof1);
    globdat.cons.CMat2 = speye(globdat.cons.nDof2);

    globdat.cons.CMat1(:,gdofs) = [];


else
    globdat.cons.nDof  = globdat.ndofs * globdat.mesh.info.nnodes;
    globdat.cons.CMat  = speye(globdat.cons.nDof); 

    globdat.cons.CMat(:,gdofs) = [];
end

% Allocate struct for info
globdat.info.errVec = cell(props.ts.nSteps,1);

end