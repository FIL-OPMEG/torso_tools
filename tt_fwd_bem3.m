function L = tt_generate_fwds_bem3(S);
% Generate the lead fields for a 3-shell bem of the chest;

if ~isfield(S,'pos'); error('please specify the source positions!'); end
if ~isfield(S,'posunits'); error('please specify the current positiosn units!');end
if ~isfield(S,'ori'); S.ori = []; end
if ~isfield(S,'sensors'); error('please specify the sensor structure!'); end
if ~isfield(S,'T'); S.T = eye(4); end
if ~isfield(S,'names'); S.names = {'blood','lungs','torso'}; end;
if ~isfield(S,'ci'); S.ci = [.62 .05 .23]; end
if ~isfield(S,'co'); S.co = [.23 .23  0 ]; end
if ~isfield(S,'isa'); S.isa = []; end


if isempty(which('hbf_BEMOperatorsPhi_LC'));
    tt_add_bem;
end

meshes = tt_load_meshes(S.T,S.names);
[~, sf] = tt_determine_mesh_units(meshes);
bmeshes = {};

figure; hold on;
for ii = 1:numel(meshes)
    bmeshes{ii}.p = meshes{ii}.vertices/sf; % conform to meters
    bmeshes{ii}.e = meshes{ii}.faces;
    plot3(bmeshes{ii}.p(:,1),bmeshes{ii}.p(:,2),bmeshes{ii}.p(:,3),'m.')
end

% get sensors into meters
S.sensors = ft_convert_units(S.sensors,'m');

% get the source space into meters
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

plot3(src.pos(:,1),src.pos(:,2),src.pos(:,3),'g*')
plot3(coils.p(:,1),coils.p(:,2),coils.p(:,3),'ko')

% conductivity for, cardiac blood, lungs, thorax
ci = S.ci;
co = S.co;


% Generate transfer matrix
D = hbf_BEMOperatorsPhi_LC(bmeshes);
if isempty(S.isa)
    Tphi_full = hbf_TM_Phi_LC(D,ci,co);
else
    fprintf('%-40s: %30s\n','Applying ISA', S.names{S.isa});
    Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,S.isa);
end
DB = hbf_BEMOperatorsB_Linear(bmeshes,coils);
TB = hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);

if isempty(S.ori)
    L = S.sensors.tra*hbf_LFM_B_LC(bmeshes,coils,TB,src.pos);
else
    L = S.sensors.tra*hbf_LFM_B_LC(bmeshes,coils,TB,S.pos,S.ori);
end