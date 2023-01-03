function props = inputData()

% initialize props as a struct
props = struct;

% Mesh file
props.meshFile = 'mesh.mphtxt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% FE model %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

props.feModel.type     = 'PhaseFieldFracture';
props.feModel.dofs     = {'dx','dy','phi'};
props.feModel.hist     = {};
props.feModel.elType   = 'tri';
props.feModel.Ipoints  = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Material Parameters %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plane_stress = false;
E            = 210e3;
nu           = 0.3;
lambda       = E*nu/((1+nu)*(1-2*nu));
mu           = E/(2*(1+nu));

if plane_stress == true
    lambda = (2*lambda*mu)/(lambda+2*mu);
end

mat        = struct;
mat.E      = E;
mat.nu     = nu;
mat.lambda = lambda;
mat.mu     = mu;
mat.K      = lambda + 2/3*mu;

mat.D      = zeros(3,3);
mat.D(1,1) = lambda + 2*mu;
mat.D(2,2) = lambda + 2*mu;
mat.D(3,3) = mu;
mat.D(1,2) = lambda;
mat.D(2,1) = lambda;

mat.gc     = 2.7; 
mat.l0     = 0.015; 
mat.Esplit = 'NoSplit'; % 'NoSplit', 'Spectral2D', 'Amor', 'Rankine'
mat.pen    = 6500;

% AT2 model
mat.p       = 2.0;
mat.a1      = 2.0;
mat.a2      = -0.5;
mat.a3      = 0.0;
mat.calpha  = 2.0;
mat.eta     = 0;

props.mat   = mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Should be struct format

cons(1).edge = {'Bottom'};
cons(1).dofs = {'dx','dy'};
cons(1).vals = {0,0};

cons(2).edge = {'Top'};
cons(2).dofs = {'dx','dy'};
cons(2).vals = {0,1};

props.cons   = cons;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Time-stepping Parameters %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts.nSteps         = 250; 
ts.initStepSize   = 1e-4;
ts.finalStepSize  = 1e-8;
ts.adaptSteps     = [56];
ts.adaptSize      = 0.02;
ts.stopRatio      = 0.05;

props.ts          = ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Solver Parameters %%%%%%¤%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlSolver.type        = 'newtonEx'; %'staggered';
nlSolver.tol         = 1e-3;
nlSolver.maxiter     = 5000;

linSolver.type       = 'direct';
linSolver.tol        = []; % Tolerance for iterative solvers
linSolver.prec       = []; % Preconditioner type for iterative solvers
linSolver.restart    = []; % Restart for iterative solvers
linSolver.maxiter    = []; % Maximum iterations for iterative solvers

props.nlSolver       = nlSolver;
props.linSolver      = linSolver;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Post-processing Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

props.postProc.printVTK    = 1;  % Prints VTK file every 'n' steps

end