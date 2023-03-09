function [meshes, names] = tt_load_meshes(T)

% Load meshes, optional argument to give affine matrix to transform;

if ~nargin
    T = eye(4);
end

meshes = {};
names = {'blood','lungs','torso'};

for ii = 1:numel(names);
    
    meshes{ii} = export(gifti(fullfile(tt_path,[names{ii} '.gii'])),...
        'patch');
    meshes{ii} = spm_mesh_transform(meshes{ii},T);
    meshes{ii}.faces = double(meshes{ii}.faces);
    
end

    