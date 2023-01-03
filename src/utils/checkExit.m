function globdat = checkExit(props,globdat)

if globdat.ts.step == props.ts.nSteps || globdat.lodi(globdat.ts.step+1,2)/...
                      max(globdat.lodi(:,2)) < props.ts.stopRatio || ...
                      globdat.ts.dt < props.ts.finalStepSize
    fprintf('\n')
    disp('*** Terminating the solution process! ***')

    globdat.active = false;
end


end