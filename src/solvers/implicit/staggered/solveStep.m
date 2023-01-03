%-------------------------------------------------------------------------%
% globdat = solveStep(props,globdat) solves a step of the nonlinear problem
%           with a staggered (alternate minimization) solver.
%
% INPUT:  props   -> struct that stores all input data
%         globdat -> global database struct
%
% OUTPUT: globdat -> globdat database struct (updated)
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   18.12.2022
%-------------------------------------------------------------------------%

function globdat = solveStep(props,globdat)

% Increment step and time
globdat.ts.step = globdat.ts.step + 1;

% Check if step-size needs to be changed
if ismember(globdat.ts.step,props.ts.adaptSteps) == 1
    globdat.ts.dt = globdat.ts.dt*props.ts.adaptSize;
end

% Update current (pseudo) time
globdat.ts.t    = globdat.ts.t + globdat.ts.dt;

% Print some information to the command window
fprintf('\n')
disp('-------------------------------------------------------')
disp(['Step ', num2str(globdat.ts.step),...
      ', t = ',num2str(globdat.ts.t), ...
      ', dt = ',num2str(globdat.ts.dt)])
disp(['Elements = ',num2str(globdat.mesh.info.nelems.(props.feModel.elType)),...
      ', DOFs = ',num2str(length(globdat.state))])
disp('-------------------------------------------------------')

% Setup error vector of size max iterations
errVec = zeros(props.nlSolver.maxiter,1);

% Set iteration count to zero and apply Dirichlet bc
iter = 0;
globdat.state = setConstraints(globdat.state,globdat.ts.t,...
                               globdat.cons.tab);

% Begin staggered iterations
while true

    % Increase iteration count
    iter = iter + 1;

    % Assemble first system of equations
    globdat  = ModelAction(props,globdat,'GET_MATRIX0_1');

    % Solve for solution increment
    dstate1   = (globdat.cons.CMat1' * globdat.K1 * ...
               globdat.cons.CMat1)\(globdat.cons.CMat1' * -globdat.fint1);
    dstate1   = globdat.cons.CMat1 * dstate1;

    % Update solution vector
    globdat.state(1:globdat.cons.nDof1) = globdat.state(1:globdat.cons.nDof1) + ...
               dstate1;

    % Assemble second system of equations
    globdat  = ModelAction(props,globdat,'GET_MATRIX0_2');

    % Solve for solution increment
    dstate2   = globdat.K2\-globdat.fint2;

    % Update solution vector
    globdat.state(globdat.cons.nDof1+1:globdat.cons.nDof1+globdat.cons.nDof2) = ...
        globdat.state(globdat.cons.nDof1+1:globdat.cons.nDof1+globdat.cons.nDof2) + ...
               dstate2;

   % Compute error in the solution and store in error vector
   err            = norm(dstate1)+norm(dstate2);
   errVec(iter,1) = err;

   % Display current iteration and error info
   disp(['Iteration ', num2str(iter), ': Error = ', num2str(err)])

   % Check for convergence
   if err < props.nlSolver.tol

       % Display message to command window            
       fprintf('\n')
       disp(['Converged in ',num2str(iter),' iteration(s).'])
       
       % Compute load and display to command window
       loads = abs(sum(globdat.fint1(globdat.cons.ldofs)));
       disp(['Load = ',num2str(loads),' Displacement = ',num2str(globdat.ts.t)])

       % Store load-displacement data
       globdat.lodi(globdat.ts.step+1,:) = [globdat.ts.t,loads];

       % Ask Model to commit
       globdat = ModelAction(props,globdat,'COMMIT');
       
       % Re-size error vector and store in globdat
       errVec(iter+1:props.nlSolver.maxiter) = [];
       globdat.errVec{globdat.ts.step,1} = errVec;

       break;

   elseif iter == props.nlSolver.maxiter

        %-------------------------------------------------------------%
        % Terminate iterations when max iter is reached or when       %
        % the error blows up                                          %
        %-------------------------------------------------------------%
        elseif iter == props.nlSolver.maxiter || err > 5.0
            
            % Message to command window
            disp(['*** WARNING: Failed to converge in ',num2str(iter),' iterations! ***'])
            disp(['*** Reducing stepsize by factor ',num2str(props.ts.adaptSize)])
            
            % Revert variables to oldStep
            globdat.state  = globdat.state0;
            
            % Reduce the stepsize for prescribed displacement method
            globdat.ts.step = globdat.ts.step - 1;
            globdat.ts.t    = globdat.ts.t - globdat.ts.t;
            globdat.ts.dt   = globdat.ts.dt*props.ts.adaptSize;
                    
            break;
   end

end

end