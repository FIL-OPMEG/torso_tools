function [meshes, names] = tt_load_meshes_grb(T,names)

% Load meshes, optional argument to give affine matrix to transform;
% optional names allows one to specify other gifti meshes in tt directory

if nargin<2,
    names=[];
end;

if nargin<1
    T=[];
end;

if isempty(T),
    T = eye(4);
end;

if isempty(names)
    names = {'blood','lungs','torso'};
end;

meshes = {};

for ii = 1:numel(names);
    
    meshes{ii} = export(gifti(fullfile(tt_path,[names{ii} '.gii'])),...
        'patch');
    meshes{ii} = spm_mesh_transform(meshes{ii},T);
    meshes{ii}.faces = double(meshes{ii}.faces);
    
end

    