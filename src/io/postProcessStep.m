function globdat = postProcessStep(dir_output,props,globdat)

if rem(globdat.ts.step,props.postProc.printVTK) == 0 || ...
           globdat.active == 0 || globdat.ts.step == 0
    
    VTKPostProcess(dir_output, props, globdat)

end

end