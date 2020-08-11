clear;clc;
%% parameters
modelname = 'model-1';
np = 300;                           % number of particles
L = 1;                              % length of bounding box
W = 1;                              % width of bounding box
H = 1;                              % height of bounding box
lc = 0.4;                           % cut-off length of connection
M = 1.0;                            % mass of particle (mean)
msd = 0.1;                          % mass of particle (stdev)
R = 1.0;                            % radius of particle (mean)
rsd = 0.1;                      	% radius of particle (stdev)
K = 1000;                           % spring stiffness (mean)
ksd = 10.0;                         % spring stiffness (stdev)

%% generate points
x = L*rand(np, 1);   %x-coordinate of a point
y = W*rand(np, 1);   %y-coordinate of a point
z = H*rand(np, 1);   %z-coordinate of a point
p = [x y z];
mass = normrnd(M, msd, [np,1]);     % mass array
radius = normrnd(R, rsd, [np,1]);	% radius array
inertia = 2/5*mass.*radius.^2;      % inertia array
%% generate connections
DT = delaunayTriangulation(p);
nt = size(DT.ConnectivityList,1);
ne = 4*nt;
edges = zeros(ne, 2);
for i = 1 : nt
    ele = DT.ConnectivityList(i,:);
    edges((i-1)*4+1,:) = ele([1 2]);
    edges((i-1)*4+2,:) = ele([2 3]);
    edges((i-1)*4+3,:) = ele([3 4]);
    edges((i-1)*4+4,:) = ele([4 1]);
end
edges = sort(edges,2);
edges = unique(edges,'rows');
edgeLength = vecnorm(p(edges(:,1),:)-p(edges(:,2),:),2,2); 
validEdges = find(edgeLength < lc);
edges = edges(validEdges,:);
edgeLength = edgeLength(validEdges);
springstiffness = normrnd(K, ksd, [size(edges,1),1]); % stiffness array
connCount = zeros(np,1);
for i = 1 : np
    connCount(i,1) = length(find(edges==i));
end
zeroConn = find(connCount==0);
if ~isempty(zeroConn)
    nzc = length(zeroConn);
    disp(['Warning: there are ', num2str(nzc), ' particles not connected with any other particles']);
end
%% generate output and visualization
% node
fid = fopen([modelname,'.geometry'],'w');
nodes = p;
[nn, ~] = size(nodes);
fprintf(fid,'%d\n', nn);
for i = 1 : nn
    fprintf(fid,'%d\t%16.8e\t%16.8e\t%16.8e\t%16.8e\t%16.8e\n', ...
        i, nodes(i,:), mass(i), inertia(i));
end
% element
elements = edges;
[ne, ~] = size(elements);
fprintf(fid,'%d\n', ne);
for i = 1 : ne
    fprintf(fid,'%d\t%d\t%d\t%16.8e\n', i, elements(i,:), springstiffness(i));
end
fclose(fid);
figure; box on; hold on; axis equal;
scatter3(x, y, z, 40, 'filled');
for i = 1 : size(edges,1)
    p1 = p(edges(i,1),:);
    p2 = p(edges(i,2),:);
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'k-', 'LineWidth', 0.5);
end

