function L = tt_generate_fwds_bem1(S);
% Generate the lead fields for a 1-shell bem of the chest;

if ~isfield(S,'pos'); error('please specify the source positions!'); end
if ~isfield(S,'posunits'); error('please specify the current positiosn units!');end
if ~isfield(S,'ori'); S.ori = []; end
if ~isfield(S,'sensors'); error('please specify the sensor structure!'); end
if ~isfield(S,'T'); S.T = eye(4); end


if isempty(which('hbf_BEMOperatorsPhi_LC'));
    tt_add_bem;
end

[meshes, names] = tt_load_meshes(S.T);
[~, sf] = tt_determine_mesh_units(meshes);
bmeshes = {};

bid = find(contains(names,'torso'));
bmeshes{1}.p = meshes{bid}.vertices/sf; % conform to m;
bmeshes{1}.e = meshes{bid}.faces;


% get sensors into meters 
S.sensors = ft_convert_units(S.sensors,'m');

% get the source space into meters
cfg                     = [];
cfg.method              = 'basedonpos';
cfg.sourcemodel.pos     = S.pos;
cfg.sourcemodel.unit     = S.posunits;
src                     = ft_prepare_sourcemodel(cfg);
src                     = ft_convert_units(src,'m');

coils = [];
coils.p = S.sensors.coilpos;
coils.n = S.sensors.coilori;

% conductivity for thorax
ci = .23;
co = 0;

% Generate transfer matrix
D = hbf_BEMOperatorsPhi_LC(bmeshes);
Tphi_full = hbf_TM_Phi_LC(D,ci,co);
DB = hbf_BEMOperatorsB_Linear(bmeshes,coils);
TB = hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);

if isempty(S.ori)
    L = S.sensors.tra*hbf_LFM_B_LC(bmeshes,coils,TB,src.pos);
else
    L = S.sensors.tra*hbf_LFM_B_LC(bmeshes,coils,TB,S.pos,S.ori);
end