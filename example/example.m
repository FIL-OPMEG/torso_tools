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
S.plot = 0;

T = tt_register_torso(S);

%% Get the units of the registered mesh

tmp = tt_load_meshes(T);
[unit, sf] = tt_determine_mesh_units(tmp);


%% Try and plot the model 

S = [];
S.subject = mesh;
S.T = T;
S.sensors = ft_convert_units(grad,unit); % <-  check the units

tt_check_registration(S);

%% Generate the source space

S = [];
S.subject = mesh;
S.T = T;
S.width = 80;
S.depth = 20:10:80;
S.resolution = 10;
S.mask = 1;

src = tt_generate_spine_grid(S);

scatter3(src.pos(:,1),src.pos(:,2),src.pos(:,3),'y.')


%% Forward modelling - Stenroos' BEM

src = ft_convert_units(src,'m'); % sources need to be in meters for BEM

S = [];
S.pos = src.pos;
S.T = T;
S.sensors = grad;

Ls = tt_fwd_bem3(S);