function [M]=tt_wrap_cord(sourcepos,rad)
% function wrap_cord(sourcepos,rad)
% sourcespos are positions of sources (1 at each depth) moving along cord
% rad is the radius of wrapper
% M is the wrapper surface mesh

[dum,vertind]=max(std(sourcepos)); %% get main vertical axis
[dum,sortind]=sort(sourcepos(:,vertind));
cp=sourcepos(sortind,:); % in ascending order

% Use icoseheda to generate cord volume
allp = [];
ico_prev = spm_mesh_sphere(0);
for ii = 1:size(cp,1)
    ico_curr = spm_mesh_sphere(0);
    T = rad * eye(4);
    T(:,end) = [sourcepos(ii,:) 1];
    ico_curr = spm_mesh_transform(ico_curr,T);
    [~, currID] = intersecting_points(ico_prev, ico_curr);
    tmp = ico_curr.vertices;
    tmp(currID,:) = [];
    allp = cat(1,allp,tmp);
    ico_prev = ico_curr;
end

% Look for the most complex boundary which satisfies a closed mesh to bound
% a fminsearch later
fprintf('Determining maxmimum valid boundary shrinkage...\n')
wraps = [0 1];
wrap_old = mean(wraps);
dWrap = Inf;
while dWrap >= 1e-2
    wrap = mean(wraps);
    k = boundary(allp(:,1),allp(:,2),allp(:,3),wrap);
    M.faces=k;
    M.vertices=allp;
    wrap_old = wrap;
    closed = is_closed(sourcepos,M);
    fprintf('%-40s: %30s\n',num2str(wrap),closed);
    switch closed
        case 'CLOSED'
            wraps(1) = wrap;
        case 'OPEN'
            wraps(2) = wrap;
    end
    dWrap = abs(wrap_old - mean(wraps));
end

if strcmpi(closed,'OPEN')
    maxShrink = mean(wraps);
    % perform one last iteration to make it a closed surface
    k = boundary(allp(:,1),allp(:,2),allp(:,3),maxShrink);
    M.faces=k;
    M.vertices=allp;
    fprintf('%-40s: %30s\n',num2str(maxShrink),...
        num2str(is_closed(sourcepos,M)));
else
    maxShrink = wrap_old;
end

options = optimset('Display', 'iter');
fprintf('Determining optimal shrinkage...')
fun = @(x) -1 * mesh_cost_fun(x,allp);
optShrink = fminbnd(fun,0,maxShrink,options);
k = boundary(allp(:,1),allp(:,2),allp(:,3),optShrink);
M.faces=k;
M.vertices=allp;
fprintf('Removing redundant vertices...\n')
M = cleanup_mesh(M);
fprintf('COMPLETE\n')

end

function res = is_closed(pos,M)
% check if every spine point is inside the mesh to determine its a closed mesh
for ii = 1:size(pos,1)
    inside(ii) = tt_is_inside(pos(ii,:),M.vertices,M.faces);
end
if all(inside)
    res = 'CLOSED';
else
    res = 'OPEN';
end
end

function [old_vert, new_vert] = intersecting_points(old,new)
% return ther vertex IDs of two meshes which intersect in space

inside = [];
for ii = 1:size(old.vertices,1)
    inside(ii) = tt_is_inside(old.vertices(ii,:),...
        new.vertices, new.faces);
end
old_vert = find(inside);

inside = [];
for ii = 1:size(new.vertices,1)
    inside(ii) = tt_is_inside(new.vertices(ii,:),...
        old.vertices, old.faces);
end
new_vert = find(inside);

end

function F = mesh_cost_fun(s,pts)
% cost function to determine the complexity of the mesh, it should
% minimize both the number of faces in the mesh and the edge distance
% smax is set to ensure that the optimized mesh is a closed mesh!
k = boundary(pts(:,1),pts(:,2),pts(:,3),s);
nfaces = size(k,1);
mD = 0;
for ii = 1:size(k,1)
    p = pts(k(ii,:),:);
    [~, D] = knnsearch(p,p,'k',2);
    mD = max(max(D(:)),mD);
end
F = -log(nfaces) -log(mD);
end

function Mclean = cleanup_mesh(M)
% remove the vertices which are not a part of a face.
usedVertices = unique(M.faces(:));
unusedVertices = setdiff(1:size(M.vertices, 1), usedVertices);
Mclean.vertices = M.vertices;
Mclean.vertices(unusedVertices,:) = [];
[~,ord] = sort(usedVertices,'ascend');
tmp = M.faces;
for ii = 1:numel(usedVertices)
    tmp(tmp==usedVertices(ii)) = ord(ii);
end
Mclean.faces = tmp;
end