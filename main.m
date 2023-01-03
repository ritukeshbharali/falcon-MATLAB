%=========================================================================%
% falcon-MATLAB is a finite element analysis program from the falcon
% numerical analyses suite. It is developed primarily for prototyping.
% The coding style is similar to the original falcon package, which is
% based on the Jem-Jive library (https://dynaflow.com/software/jive/).
% 
% Author   : Ritukesh Bharali (Chalmers University of Technology)
%
%=========================================================================%

function globdat = main(problem)

% User prompt if no problem is provided
if nargin ~= 1
    problem = input('Enter problem to run: ','s');
end

% Clear the workspace
close all; clc
format long;
warning('off')


%=========================================================================%
% ADD PATH TO RELEVANT DIRECTORIES
%=========================================================================%

disp(' - Adding path to directories')

% Remove path
rmpath(genpath('./src/fem'))
rmpath(genpath('./src/solvers'));
rmpath(genpath('./inputData'))

% Path to source folders
addpath(genpath('./ext'));
addpath(genpath('./src/constraints'));
addpath(genpath('./src/io'));
addpath(genpath('./src/solvers/linear'));
addpath(genpath('./src/utils'));

% Path to input files
addpath(fullfile('./inputData/',problem));

% Create an output directory
dir_output = sprintf('./output/%s', problem);
if ~exist(dir_output, 'dir')
    mkdir(dir_output)
end
addpath(dir_output)


%=========================================================================%
% INITIALIZATION
%=========================================================================%

disp(' - Reading user input and setting up runtime parameters')
props = inputData;

% Add path to the FE Model and solver
addpath(fullfile('./src/fem/',props.feModel.type))
addpath(fullfile('./src/solvers/implicit/',props.nlSolver.type))

% Initialise the global database
globdat        = initProblem(props);

% Create PVD file
globdat = postProcessStep(dir_output,props,globdat);


%=========================================================================%
% RUN
%=========================================================================%

tic

% Startup time-stepping
globdat.active = true;

while globdat.active
           
    % Solve a step
    globdat = solveStep(props,globdat);
    
    % Post-process
    globdat = postProcessStep(dir_output,props,globdat);

    % Check for exit
    globdat = checkExit(props,globdat);
    
end

%=========================================================================%
% SHUTDOWN
%=========================================================================%

% Shutdown problem
globdat     = shutdownProblem(dir_output,props,globdat);

toc

end