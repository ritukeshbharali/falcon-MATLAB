%-------------------------------------------------------------------------%
% globdat = solveStep(props,globdat) solves a step of the nonlinear problem
%           with a monolithic Newton-Raphson solver.
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
% DATE:   20.12.2022
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
    globdat  = ModelAction(props,globdat,'GET_MATRIX0');

    % Solve for solution increment
    dstate   = (globdat.cons.CMat' * globdat.K * ...
               globdat.cons.CMat)\(globdat.cons.CMat' * -globdat.fint);
    dstate   = globdat.cons.CMat * dstate;

    % Update solution vector
    globdat.state = globdat.state + dstate;

    % Compute error in the residual and store in error vector
    if iter == 1
        err0 = norm((globdat.cons.CMat' * globdat.fint));
        err  = 1;
    else
        err = norm(globdat.cons.CMat'*globdat.fint)/max(1,err0);
    end
    errVec(iter,1) = err;

    % Display current iteration and error info
    disp(['Iteration ', num2str(iter), ': Error = ', num2str(err)])

    % Check for convergence
    if err < props.nlSolver.tol

       % Display message to command window            
       fprintf('\n')
       disp(['Converged in ',num2str(iter),' iteration(s).'])
       
       % Compute load and display to command window
       loads = abs(sum(globdat.fint(globdat.cons.ldofs)));
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
        elseif iter == props.nlSolver.maxiter || err > 1e4
            
            % Message to command window
            disp(['*** WARNING: Failed to converge in ',num2str(iter),' iterations! ***'])
            disp(['*** Reducing stepsize by factor ',num2str(props.ts.adaptSize)])
            
            % Revert variables to oldStep
            globdat.state  = globdat.state0;
            
            % Reduce the stepsize for prescribed displacement method
            globdat.ts.step = globdat.ts.step - 1;
            globdat.ts.t    = globdat.ts.t - globdat.ts.t;
            globdat.ts.dt   = globdat.ts.dt*props.ts.adaptSize;

            % Ask Model to revert
            globdat = ModelAction(props,globdat,'REVERT');
                    
            break;
    end

end

end