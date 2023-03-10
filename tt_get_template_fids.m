function fids = tt_get_template_fids(T)

if nargin == 0
    T = eye(4)
end

[meshes, names] = tt_load_meshes(T);

id = find(contains(names,'torso'));

% template fiducuals of the torso
% - left shoulder   (pt 901)
% - right shoulder  (pt 1953)
% - low spine       (pt 1147)

fids = meshes{id}.vertices([901 1953 1147],:);


