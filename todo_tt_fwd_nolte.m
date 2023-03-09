function L = tt_fwd_nolte(S);
% % Generate the lead fields for a nolte corrected sphere of the chest;
% 
% if ~isfield(S,'pos'); error('please specify the source positions!'); end
% if ~isfield(S,'ori'); S.ori = []; end
% if ~isfield(S,'sensors'); error('please specify the sensor structure!'); end
% if ~isfield(S,'T'); S.T = eye(4); end
% 
% 
% [meshes, names] = tt_load_meshes(S.T);
% bmeshes = {};
% 
% bid = find(contains(names,'thorax'));
% bmeshes{1}.p = meshes{bid}.vertices;
% bmeshes{1}.e = meshes{bid}.faces;
% 
% % get sensors into meters 
% S.sensors = ft_convert_units(S.sensors,'m');
% 
% % use the fitted thorax
% bnd = [];
% bnd.tri         = bmeshes{id}.e;
% bnd.pnt         = bmeshes{id}.p;
% bnd.unit        = 'm';
% 
% cfg                     = [];
% cfg.method              = 'singleshell';
% vol                     = ft_prepare_headmodel(cfg,bnd);
% 
% % source space only using the points inside the bem;
% cfg                     = [];
% cfg.method              = 'basedonpos';
% cfg.sourcemodel.pos     = fwd.pos(inid,:)/1000;
% cfg.sourcemodel.inside  = ones(length(cfg.sourcemodel.pos),1);
% src                     = ft_prepare_sourcemodel(cfg);
% 
% % calculate forwards
% cfg                     = [];
% cfg.sourcemodel         = src;
% cfg.headmodel           = vol;
% cfg.grad                = grad;
% % cfg.reducerank          = 'yes';
% 
% fwd_ss_thorax = ft_prepare_leadfield(cfg);
% 
% if isempty(S.ori)
%     L = S.sensors.tra*hbf_LFM_B_LC(bmeshes,coils,TB,src.pos);
% else
%     L = S.sensors.tra*hbf_LFM_B_LC(bmeshes,coils,TB,S.pos,S.ori);
% end