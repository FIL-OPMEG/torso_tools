function tt_check_registration(S)

if ~isfield(S,'subject'); error('please provide a subject mesh'); end
if ~isfield(S,'T'); warning('no transformation matrix found'); S.T = eye(4); end
if ~isfield(S,'sensors') S.sensors = []; end

figure
ft_plot_mesh(S.subject,'facecolor','none','edgecolor','k',...
    'clipping','off','edgealpha',0.5);

hold on

colours = {'r','g','b'};
alpha = [0.3 0.3 0.3];

meshes = tt_load_meshes(S.T);
unit = tt_determine_mesh_units(meshes);

for ii = 1:numel(meshes);
    
    tmp = [];
    tmp.pnt = meshes{ii}.vertices;
    tmp.tri = meshes{ii}.faces;
    tmp.unit = unit;
    
    ft_plot_mesh(tmp,'facecolor',colours{ii},'edgecolor','none','facealpha',alpha(ii))
    
    
end

if ~isempty(S.sensors)
    ft_plot_sens(ft_convert_units(S.sensors,unit))
end