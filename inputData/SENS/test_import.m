clear all
close all
clc

%----------------------------------------------------%
disp('Trying to import a COMSOL mphtxt mesh ...')

file_name = 'mesh.mphtxt';
mesh = readCOMSOLmphtxt(file_name);

nodes = mesh.nodes;
elems = mesh.elems;
idx   = mesh.idx;
info  = mesh.info;

% disp('Trying to import a COMSOL mphtxt mesh ...')
% 
% % Extract top edge elements
% top_el    = find(ismember(idx.edg, idx.Top) == 1);
% top_nodes = elems.edg(top_el,2:end);
% top_nodes = unique(top_nodes);
% 
% % Extract bottom edge elements
% bot_el    = find(ismember(idx.edg, idx.Bottom) == 1);
% bot_nodes = elems.edg(bot_el,2:end);
% bot_nodes = unique(bot_nodes);


% ------ Plotting (works only for Tri3) --------- %

nelems   = info.nelems.tri;
nelnodes = length(elems.tri(1,:))-1;

X = zeros(nelnodes,nelems);
Y = zeros(nelnodes,nelems);

for i = 1:nelems
    node   = elems.tri(i,2:nelnodes+1);
    X(:,i) = nodes(node,2);
    Y(:,i) = nodes(node,3);
end

figure
fill(X,Y,'w')
hold on
axis square

%pause(1)

%close all

disp('Success! ')



