function headmodel = ft_headmodel_bemhbf(mesh, varargin)

% FT_HEADMODEL_BEMHBF generates the BEM solution for M/EEG conductive
% modelling 

% add toolbox
if isempty(which('hbf_BEMOperatorsPhi_LC'))
    tt_add_bem;
end

% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');
isolated        = ft_getopt(varargin, 'isolatedsource', []);

% conductivity values inside and outside of tissue types
ci = conductivity(1,:);
co = conductivity(2,:);

% convert meshes into hbf_bem compatible versions
bmeshes = [];
for ii = 1:numel(mesh)
    bmeshes{ii}.p = mesh(ii).pos;
    bmeshes{ii}.e = mesh(ii).tri;
end

% Generate transfer matrix
D = hbf_BEMOperatorsPhi_LC(bmeshes);
if isempty(isolated)
    Tphi_full = hbf_TM_Phi_LC(D,ci,co);
else
    fprintf('%-40s: %30s\n','Applying ISA to compartment',...
        num2str(isolated));
    Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,isolated);
end

% pack everything up
headmodel           = [];
headmodel.bmeshes   = bmeshes;
headmodel.cond      = conductivity;
headmodel.mat       = Tphi_full;
headmodel.type      = 'bem_hbf';
headmodel.unit      = mesh(1).unit;
