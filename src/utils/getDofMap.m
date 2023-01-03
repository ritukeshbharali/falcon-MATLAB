function dof = getDofMap(nodes,dim)

n    = length(nodes)*dim;
dof  = zeros(n,1);

switch dim

    case 2
        dof(1:2:n) = 2*nodes-1;
        dof(2:2:n) = 2*nodes;

    case 3
        dof(1:3:n) = 3*nodes-2;
        dof(2:3:n) = 3*nodes-1;
        dof(3:3:n) = 3*nodes;
    otherwise
        error('Choose 2 for 2D!')
end

end

