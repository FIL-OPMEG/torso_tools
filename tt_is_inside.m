function inside = tt_is_inside(pos,mesh_pos,mesh_tri)

% check the mesh faces and positions are all doubles otherwise mex file
% will crash
pos = double(pos);
mesh_pos = double(mesh_pos);
mesh_tri = double(mesh_tri);

sa = abs(sum(solid_angle(mesh_pos-pos,mesh_tri)));

if abs(sa - 4*pi) < 1e-6
    inside = 1;
else
    inside = 0;
end