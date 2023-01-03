function globdat = shutdownProblem(dir_output,props,globdat)

% Re-size lodi array and export to txt file
globdat.lodi(globdat.ts.step+2:end,:) = [];
writematrix(globdat.lodi,fullfile(dir_output,'lodi.txt'))  

end

