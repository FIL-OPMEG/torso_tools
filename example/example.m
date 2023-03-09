clearvars
close all
clc

% Should run with the FIL copy of SPM

% addpath D:\Documents\MATLAB\torso_tools\

load(fullfile(tt_path,'example','sensors.mat'));

%% Register the thorax to subject model

fname = fullfile(tt_path,'example','seated_body_registered.stl');

subject = ft_read_headshape(fname);
p = [];
p.vertices = subject.pos;
p.faces = subject.tri;
p2 = reducepatch(p,0.6);
subject.pos = p2.vertices;
subject.tri = p2.faces;

mesh = [];
mesh.vertices = subject.pos;
mesh.faces = subject.tri;


% spm_mesh_render(p2);

% get the fiducuals of the torso
% - left shoulder   (pt ?)
% - right shoulder  (pt ?)
% - low spine       (pt ?)
sub_fids = [1072 -614 161
    618 -569 145
    864 -466 -330];
% sub_fids = [1072 -614 161
%     618 -569 145
%     864 -480 -164];



S = [];
S.subject = mesh; 
S.fiducials = sub_fids;
S.plot = 1;

T = tt_register_torso(S);


%% Try and plot the model 

S = [];
S.subject = mesh;
S.T = T;
S.sensors = ft_convert_units(grad,'mm'); % <-  check the units

tt_check_registration(S);

%% Determine which sources are inside the BEM to compare + plot

[bmeshes_reg, names] = tt_load_meshes(T);

id = find(contains(names,'torso'));


src = ft_convert_units(src,'mm');
grad = ft_convert_units(grad,'mm');

for ii = 1:length(src.pos)

    
    tmp = src.pos(ii,:); 
    inside(ii) = tt_is_inside(tmp,...
        bmeshes_reg{id}.vertices,bmeshes_reg{id}.faces);
    
    
end

sum(inside)
inid = find(inside);

scatter3(src.pos(inid,1),src.pos(inid,2),src.pos(inid,3),'y*')


%% Forward modelling - Stenroos' BEM


S = [];
S.pos = src.pos(inid,:);
S.T = T;
S.sensors = grad;

Ls = tt_fwd_bem3(S);