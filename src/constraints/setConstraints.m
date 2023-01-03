%-------------------------------------------------------------------------%
% state = setConstraints(state,t,cons) applies time-dependent Dirichlet 
% Boundary Condition to the solution vector 'state'. 
%
% INPUT:  state   -> Solution vector
%         t       -> Current time
%         cons    -> Matrix containing [DirichletDOFS; Factor]
%                    Factor = 0 : homogeneous Dirichlet condition
%                    Factor = 1 : ihhomogeneous Dirichlet condition
%
%
% OUTPUT: state   -> Solution vector with Dirichlet values applied
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   17.12.2022
%-------------------------------------------------------------------------%

function state = setConstraints(state,t,cons)

dofs           = cons(:,1);
fact           = cons(:,2);
state(dofs,1)  = t * fact;

end