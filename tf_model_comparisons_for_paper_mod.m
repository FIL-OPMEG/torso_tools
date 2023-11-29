clearvars
close all
clc

% if isempty(which('ft_defaults'))
%     go_fieldTrip
% end

% Comparison between Matti's BEM model of the chest and the Nolte's Single
% Shell, which we've been using.

torsotoolsdir='D:\torso_tools\'
addpath(torsotoolsdir)

sens_type = 'dense';

CERVICAL=1; %% cervical part of cord
FRONT=1; %% front of torso . These options only useful for sens_type='dense'
shift = 0;

files.results = fullfile('.','figures_real_geometries',['sens_' sens_type]);
if ~exist(files.results,'dir');
    mkdir(files.results)
end

%% Register the thorax to subject model

subject = ft_read_headshape('D:\torso_tools\example\seated_body_registered.stl');
p = [];
p.vertices = subject.pos;
p.faces = subject.tri;
p2 = reducepatch(p,0.6);
subject.pos = p2.vertices;
subject.tri = p2.faces;
subject.unit = 'mm';


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


S = [];
S.subject = mesh; % must be in units of m
S.fiducials = sub_fids;
S.plot = 0;
T = tt_register_thorax(S);

%% Load sensors



switch sens_type
    case 'real'
        load([torsotoolsdir '\example\sensors.mat']);
    case 'dense'
        S = [];
        S.T = T;
        S.subject = mesh;
        S.resolution=30; %% sensor spacing
        S.frontflag=FRONT; % look at front of torso
        if CERVICAL, %% only cervical portion of cord
            S.zlim=[100 Inf];
        end;
        grad = tt_generate_sensor_array(S);
    otherwise
        error('not an option')
end
L = [];

%% Shift sensors


id = find(contains(grad.label,'-Y'));

vec = mean(grad.coilori(id,:));

% vec = [-0.0938    0.9955    0.0113];
vec = vec./norm(vec);

grad.coilpos = grad.coilpos - shift*vec;
grad.chanpos = grad.chanpos - shift*vec;

%% Try and plot the model

S = [];
S.subject = mesh;
S.T = T;
S.sensors = grad;

tt_check_registration(S);

%% Generate source space grid + plot

S = [];
S.subject = mesh;
S.T = T;
S.width = 1;
S.depth = 50;
S.resolution = 30;
S.mask = 1;
if CERVICAL, %% only cervical portion of cord
    S.zlim=[100 +Inf];
end;

sources = tt_generate_spine_grid(S);
myrad=8; %% radius of enclosing cylinder in mm

Mw=tt_wrap_cord(sources.pos,myrad); %% create mesh around cord of radius myrad
shiftpos=Mw.vertices-repmat([0 30 0],length(Mw.vertices),1).*0;
%Mw.vertices=shiftpos+randn(length(shiftpos),3);

%
scatter3(sources.pos(:,1),sources.pos(:,2),sources.pos(:,3),'y.')

%% For paper, plot sensor layouts

figure(1000);clf
ft_plot_mesh(subject,'facecolor',[0.6 0.6 0.6],'edgecolor','none','facealpha',0.5)
hold on
scatter3(sources.pos(:,1),sources.pos(:,2),sources.pos(:,3),'r.')
ft_plot_sens(ft_convert_units(grad,subject.unit),'coilshape','sphere','facecolor','k','coilsize',15)

trisurf(Mw.faces,Mw.vertices(:,1),Mw.vertices(:,2),Mw.vertices(:,3))

set(gcf,'color','w')
view([82 0])
saveas(gcf,sprintf('source_sensors_%s_side.pdf',sens_type))
view([172 0])
saveas(gcf,sprintf('source_sensors_%s_real.pdf',sens_type))

%% Cylinder source model


% load D:\thorax_model\source_sensors.mat
%
% [bmeshes_reg, names] = tt_load_meshes(T);
% unit = tt_determine_mesh_units(bmeshes_reg);
%
% sources = [];
% sources.tri = spine.e;
% sources.pos = spine.p;
% sources.nrm = spine.nn;
% sources = ft_transform_geometry(T,sources);
% sources.units = unit;
%
% scatter3(sources.pos(:,1),sources.pos(:,2),sources.pos(:,3),'y.')

%% Determine which sources are inside the BEM to compare + plot
%
% [bmeshes_reg, names] = tt_load_meshes(T);
% [unit, sf] = tt_determine_mesh_units(bmeshes_reg);
%
% id = find(contains(names,'thorax'));
%
% sources = ft_convert_units(sources,unit);
%
% for ii = 1:length(sources.pos)
%
%
%     tmp = sources.pos(ii,:);
%     inside(ii) = tt_is_inside(tmp,...
%         bmeshes_reg{id}.vertices,bmeshes_reg{id}.faces);
%
%
% end
%
%
% inid = find(inside);


%% Determine planar distace, keep everything more than 6mm from a boundary

% [unit,sf] = tt_determine_mesh_units(bmeshes_reg)
%
% for jj = 1:3
%
%     m = [];
%     m.faces = bmeshes_reg{jj}.faces;
%     m.vertices = bmeshes_reg{jj}.vertices;
%
%     [~,norms] = spm_mesh_normals(m,1);
%
%     idx = knnsearch(m.vertices,sources.pos);
%
%
%     for ii = 1:length(sources.pos);
%
%         X = sources.pos(ii,:);
%         % Find the nearest vertex
%         p0 = m.vertices(idx(ii),:);
%
%         % Find all faces which contain this vertex
%         fid = find(sum(m.faces==idx(ii),2));
%         nrms = norms(fid,:);
%
%         D = (X-p0)*nrms';
%         Dplane(ii,jj) = min(abs(D));
%
%
%     end
%
% end
%
% inid = find(sum(Dplane<=0.006/sf,2)==0);
%
% % inid_final = intersect(inid,farid);
%
% scatter3(sources.pos(inid,1),sources.pos(inid,2),...
%     sources.pos(inid,3),'y*')
% drawnow
%


%% Prepare the common assets across all forward models now

% prep source space
cfg                     = [];
cfg.method              = 'basedonpos';
cfg.sourcemodel.pos     = sources.pos;
cfg.sourcemodel.inside  = ones(length(cfg.sourcemodel.pos),1);
cfg.sourcemodel.unit    = sources.unit;
src                     = ft_prepare_sourcemodel(cfg);
src                     = ft_convert_units(src,'m');
grad                    = ft_convert_units(grad,'m');

%% Forward modelling - Infinite Medium

L = [];

cfg                     = [];
cfg.method              = 'infinite';
vol                     = ft_prepare_headmodel(cfg);
vol.type                = 'infinite_currentdipole';
vol.unit                = 'm';

% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat(fwd_tmp.leadfield);

% figure(numel(L));clf
% ft_plot_mesh(ft_convert_units(subject,'m'),'facecolor','none', ...
%     'edgecolor',[0 2 74]/255,'edgealpha',1)
% ft_plot_headmodel(vol,'facecolor','b','facealpha',0.5,'edgecolor','none')
% set(gcf,'Color',[124 124 255]/255)
% view([90 0])

%% Forward modelling - BIG SPHERE (fitted to body scan)

% use the subject model
bnd = [];
bnd.tri         = subject.tri;
bnd.pos         = subject.pos;
bnd.unit        = subject.unit;
bnd             = ft_convert_units(bnd,'m');

cfg                     = [];
cfg.method              = 'singlesphere';
vol                     = ft_prepare_headmodel(cfg,bnd);

% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat(fwd_tmp.leadfield);

% figure(numel(L));clf
% ft_plot_mesh(ft_convert_units(subject,'m'),'facecolor','none', ...
%     'edgecolor',[0 0 0],'edgealpha',0.5)
% ft_plot_headmodel(vol,'facecolor','b','facealpha',0.3,'edgecolor','none')
% set(gcf,'Color','w')
% view([90 0])
%% Forward modelling - LITTLE SPHERE (fitted to fitted throax)

[bmeshes_reg,names] = tt_load_meshes(T);
id = find(contains(names,'torso'));

bnd = [];
bnd.tri         = bmeshes_reg{id}.faces;
bnd.pos         = bmeshes_reg{id}.vertices;
bnd.unit        = tt_determine_mesh_units(bmeshes_reg);
bnd             = ft_convert_units(bnd,'m');

cfg                     = [];
cfg.method              = 'singlesphere';
vol                     = ft_prepare_headmodel(cfg,bnd);

% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = grad;
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat(fwd_tmp.leadfield);

% figure(numel(L));clf
% ft_plot_mesh(ft_convert_units(subject,'m'),'facecolor','none', ...
%     'edgecolor',[0 0 0],'edgealpha',0.5)
% ft_plot_headmodel(vol,'facecolor','b','facealpha',0.3,'edgecolor','none')
% set(gcf,'Color','w')
% view([90 0])
%% Forward modelling - Nolte on Thorax

id = find(contains(names,'torso'));

bnd = [];
bnd.tri         = bmeshes_reg{id}.faces;
bnd.pos         = bmeshes_reg{id}.vertices;
bnd.unit        = tt_determine_mesh_units(bmeshes_reg);
bnd             = ft_convert_units(bnd,'m');

cfg                     = [];
cfg.method              = 'singleshell';
vol                     = ft_prepare_headmodel(cfg,bnd);

% calculate forwards
cfg                     = [];
cfg.sourcemodel         = src;
cfg.headmodel           = vol;
cfg.grad                = ft_convert_units(grad,'m');
cfg.reducerank          = 'no';

fwd_tmp = ft_prepare_leadfield(cfg);

L{end+1} = cell2mat({fwd_tmp.leadfield{:}});

%% Forward modelling - 1/3 shell BEM

S = [];
S.pos = src.pos;
S.posunits = src.unit;
S.T = T;
S.sensors = grad;

L{end+1} = tt_fwd_bem1(S);
L{end+1} = tt_fwd_bem3(S);



MwT = spm_mesh_transform(Mw,pinv(S.T)); %% inverse transform (gets transformed back)
save(gifti(MwT),[torsotoolsdir 'whitematter.gii']) % put in mesh directory


S.names = {'blood','lungs','torso','whitematter'}
S.ci = [.62 .05 .23 0.023 ];
S.co = [.23 .23  0 0.23];
%S.ci = [ 0.1 .23];
%S.co = [  0.23 0];
L{end+1} = tt_fwd_bem3(S);

%% Do the results not assuming spltting the dipoles up

models = {'Infinite','Large Sphere','Small Sphere','Single Shell',...
    'BEM 1C','BEM 3C','BEM 4C'};
clear cc
for ii = 1:numel(L)
    for jj = 1:numel(L)

        La = L{ii};
        %         nLa = vnorm(La,1);
        %         nLa(nLa == 0) = mean(nLa(nLa ~= 0));
        Lb = L{jj};
        %         nLb(nLb == 0) = mean(nLb(nLb ~= 0));

        re = vnorm(Lb-La,1)./vnorm(La,1);
        ra = log10(vnorm(Lb,1)./vnorm(La,1));
        %         ra(isnan(ra)) = 0;
        La_mc = La - mean(La);
        Lb_mc = Lb - mean(Lb);

        for kk = 1:size(Lb,2);
            cc(kk) = abs(corr(La(:,kk),Lb(:,kk)));
            beta(kk) = La(:,kk)'*pinv(Lb(:,kk)');
        end
        %         cc(isnan(cc)) = [];

        racell{ii,jj} = median(ra,'omitnan');
        recell{ii,jj} = median(re,'omitnan');
        cccell{ii,jj} = median(cc,'omitnan');
        bcell{ii,jj} = median(abs(beta),'omitnan');

    end
end

ccmat = cell2mat(cccell);

% save(sprintf('ccmat_%02d_mm',round(1000*shift)),'ccmat');

figure(100);clf
mat = cell2mat(recell);
imagesc(mat)
axis equal
axis off
colorbar
%colormap(brewermap(100,'RdPu'))
% xline([3.5 6.5 9.5],'--')
% yline([3.5 6.5 9.5],'--')
% xlim([0.5 12.5])
caxis([0 2.3]);

for ii = 1:numel(mat)
    [x,y] = ind2sub(size(mat),ii);
    if x ~= y

        if abs(mat(ii)) < 1.3;
            col = 'k';
        else
            col = 'w';
        end
        text(y,x,sprintf('%.2f',mat(ii)),...
            'horizontalalignment','center','color',col,...
            'fontsize',14)
    elseif mat(ii) > 0
        text(y,x,sprintf('%1.e',mat(ii)),...
            'horizontalalignment','center','color','w',...
            'fontsize',12)
    end


end
set(gcf,'color','w')
% xticklabels(models);
set(gcf,'position',[   481   411   879   527]);

fname = fullfile(files.results,...
    're_mat_all.eps');
saveas(gcf,fname,'epsc');

figure(400);clf
mat = cell2mat(cccell);
imagesc(mat)
axis equal
axis off
colorbar
%colormap(brewermap(100,'Reds'))
% xline([3.5 6.5 9.5],'--')
% yline([3.5 6.5 9.5],'--')
% xlim([0.5 12.5])
caxis([0 1])
for ii = 1:numel(mat)

    [x,y] = ind2sub(size(mat),ii);
    if x ~= y
        if abs(mat(ii)) > 0.5;
            col = 'w';
        else
            col = 'k';
        end
        text(y,x,sprintf('%.2f',mat(ii)),...
            'horizontalalignment','center','color',col,...
            'fontsize',14)

    end


end

set(gcf,'color','w')
set(gcf,'position',[   481   411   879   527]);

fname = fullfile(files.results,...
    'cc_mat_all.eps');
saveas(gcf,fname,'epsc');

%% Denrograms

models = {'Infinite','Large Sphere','Small Sphere','Single Shell',...
    'BEM 1C','BEM 3C','BEM4C'};

% Lmat = [];
% for ii = 1:8
% Lmat = [Lmat; L{ii}(:)'];
% end

mat = cell2mat(recell);
mat= 0.5*(mat + mat');

mask = find(tril(ones(numel(L)) - eye(numel(L))));

D = mat(mask)';

tree = linkage(D);
leafOrder = optimalleaforder(tree,D);
figure(50);clf
H = dendrogram(tree,'reorder',leafOrder,'ColorThreshold','default');
set(H,'linewidth',2)
set(gcf,'color','w')
grid on
xticklabels({models{leafOrder}});
xtickangle(30);

ylabel('median(Relative Error)');
axis square

fname = fullfile(files.results,...
    're_dendrogram.eps');
saveas(gcf,fname,'epsc');

mat = cell2mat(cccell);
mat= 0.5*(mat + mat');

mask = find(tril(ones(numel(L)) - eye(numel(L))));

D = 1 - mat(mask)';

tree = linkage(D);
% leafOrder = optimalleaforder(tree,D);
figure(51);clf
H = dendrogram(tree,'reorder',leafOrder,'ColorThreshold','default');
set(H,'linewidth',2)
set(gcf,'color','w')
grid on
xticklabels({models{leafOrder}});
xtickangle(30);
% set(gca,'fontname','Atkinson Hyperlegible')
ylabel('1 - median(Correlation)');
axis square


fname = fullfile(files.results,...
    'cc_dendrogram.eps');
saveas(gcf,fname,'epsc');

% end

% error('end')

%% Measure eigenspectrum

for ii = 1:numel(L)
    Ltmp = L{ii};
    tic;S = svd(Ltmp,'econ');toc
    Scum(ii,:) = cumsum(S./sum(S));
    n90(ii) = knnsearch( Scum(ii,:)',0.9);
end

cmap = [0.12156862745098039, 0.4666666666666667,  0.7058823529411765
    1.0,                 0.4980392156862745,  0.0549019607843137
    0.17254901960784313, 0.6274509803921569,  0.1725490196078431
    0.8392156862745098,  0.15294117647058825, 0.1568627450980392
    0.5803921568627451,  0.403921568627451,   0.7411764705882353
    0.5490196078431373,  0.33725490196078434, 0.2941176470588235
    0.8901960784313725,  0.4666666666666667,  0.7607843137254902
    0.4980392156862745,  0.4980392156862745,  0.4980392156862745
    0.7372549019607844,  0.7411764705882353,  0.1333333333333333
    0.09019607843137255, 0.7450980392156863,  0.8117647058823529 ];

figure;clf;
leg = [];
for ii = 1:numel(L)
    plot(Scum(ii,:),'linewidth',2,'color',cmap(ii,:))
    hold on;
    leg{end+1} = sprintf('%s = %d',models{ii},n90(ii))
end
grid on
legend(leg,'location','se')
set(gcf,'color','w')
axis square
xlabel('No. Components')
ylabel('Variance Explained')

switch sens_type
    case 'real'
        xlim([1 100])
    case 'dense'
        xlim([1 250])
end

fname = fullfile(files.results,...
    'eigenspectrum.eps');
saveas(gcf,fname,'epsc');
%% Measure the angle between the dipole orientation and the sphere radius

large_sphere_o = [0.6626 -1.6975 0.1324];   %m
small_sphere_o = [0.8502 -0.5414 -0.0634];  %m

large_source_rad = sources.pos - 1000*large_sphere_o;
small_source_rad = sources.pos - 1000*small_sphere_o;

% normalise
large_source_rad = large_source_rad ./ vnorm(large_source_rad,2);
small_source_rad = small_source_rad ./ vnorm(small_source_rad,2);

large_source_rad = repmat(large_source_rad,1,3)';
large_source_rad = reshape(large_source_rad,3,[])';
small_source_rad = repmat(small_source_rad,1,3)';
small_source_rad = reshape(small_source_rad,3,[])';

nrm = repmat(eye(3),length(sources.pos),1);

% get angle of incidence
ang_large = abs(acosd(dot(large_source_rad',nrm')));
ang_small = abs(acosd(dot(small_source_rad',nrm')));

ang_large(ang_large > 90) = 180 - ang_large(ang_large > 90);
ang_small(ang_small > 90) = 180 - ang_small(ang_small > 90);

%% box/whisker to show relationship between dipole angle and metrics

models1 = {'Infinite','Large Sphere','Small Sphere','Single Shell',...
    'BEM 1C','BEM 3C','BEM 4C'};

% models2 = {'Large Sphere','Small Sphere','Nolte','BEM 3C'};
models2 = {'Infinite','Large Sphere','Small Sphere','Single Shell',...
    'BEM 1C','BEM 3C','BEM 4C'};

id = match_str(models1,models2)

for ii = 1:numel(id)
    L2{ii} = L{id(ii)};
end
clear cc
for ii = 1:(numel(id)-1)

    Lb = L2{end};
    La = L2{ii};

    re = vnorm(Lb-La,1)./vnorm(La,1);

    for kk = 1:size(Lb,2);
        cc(kk) = corr(La(:,kk),Lb(:,kk));
        beta(kk) = La(:,kk)'*pinv(Lb(:,kk)');
    end
    %         cc(isnan(cc)) = [];

    re2bem3{ii} = re';
    cc2bem3{ii} = cc';
    % max(cc)
    % min(cc)

end
names = [];
[names{1:(3*length(sources.pos))}] = deal('Infinite');
[names{end+(1:(3*length(sources.pos)))}] = deal('Large Sphere');
[names{end+(1:(3*length(sources.pos)))}] = deal('Small Sphere');
[names{end+(1:(3*length(sources.pos)))}] = deal('Single Shell');
[names{end+(1:(3*length(sources.pos)))}] = deal('BEM 1C');


% orient = [ang_large ang_small ang_small];

re = cell2mat(re2bem3);
cc = 1 - abs(cell2mat(cc2bem3));

% [bins, edges] = discretize(orient,6);

%
% figure(150);clf
%
% g = gramm('x',names','y',re(:),'color',names');
% % g.geom_jitter('width',0.4,'edgewidth',0.01,'edgecolor','none','dodge',1,'alpha',0.05) %set edgecolor to 'none' to remove any outlines of points
% % g.stat_violin('normalization','width','width',0.5,'half',1,'dodge',2)
% g.stat_boxplot('width',2,'alpha',1,'linewidth',1.5,'drawoutlier',0)
%
% g.set_order_options('x',0,'color',0)
% g.set_color_options('map','d3_10')
% g.set_names('x','Model','y','Relative Error to 3C BEM')
% % g.axe_property('Xticklabel',{'0-15','15-30','30-45','45-60','60-75','75-90'},...
%     % 'plotboxaspectratio',[1.8 1 1]);
% g.axe_property("YLim",[0 2.5],'ygrid','on')
% g.set_text_options('base_size',12);
% g.no_legend()
%
% g.draw()
% set(gcf,'position',[726.6000  333.8000  644.8000  545.6000])
%
% fname = fullfile(files.results,...
%     're_boxwhikser.eps');
% saveas(gcf,fname,'epsc');
%
% figure(250);clf
%
%
%
% clear g
% g = gramm('x',names','y',cc(:),'color',names');
% % g.geom_jitter('width',0.4,'edgewidth',0.01,'edgecolor','none','dodge',1,'alpha',0.05) %set edgecolor to 'none' to remove any outlines of points
% % g.stat_violin('normalization','width','width',0.5,'half',1,'dodge',2)
% g.stat_boxplot('width',2,'alpha',1,'linewidth',1.5,'drawoutlier',0)
%
% g.set_order_options('x',0,'color',0)
% g.set_color_options('map','d3_10')
% g.set_names('x','Model','y','1 - abs(Correlation to 3C BEM)')
% % g.axe_property('Xticklabel',{'0-15','15-30','30-45','45-60','60-75','75-90'},...
%     % 'plotboxaspectratio',[1.8 1 1]);
% g.axe_property("YLim",[0 1],'ygrid','on')
% g.set_text_options('base_size',12);
% g.no_legend()
%
% g.draw()
% set(gcf,'position',[726.6000  333.8000  644.8000  545.6000])
%
% fname = fullfile(files.results,...
%     'cc_boxwhikser.eps');
% saveas(gcf,fname,'epsc');



%% compare field maps visually
%% Try and make a topolot
sub = [];
sub.pos = mesh.vertices;
sub.tri = mesh.faces;
sub.coordsys = 'neuromag';
sub.unit = 'mm'
grad.coordsys = 'neuromag';
grad = ft_convert_units(grad,'mm')


cfg = [];
cfg.grad = grad;
cfg.projection = 'orthographic';
% cfg.headshape = sub;
cfg.viewpoint = 'anterior';
lay = ft_prepare_layout(cfg);


[dum,srcind]=sort(dot(src.pos'-median(src.pos)',src.pos'-median(src.pos)'))

figure;
plot3(src.pos(:,1),src.pos(:,2),src.pos(:,3),'.c');
hold on;
plot3(src.pos(srcind(1:3),1),src.pos(srcind(1:3),2),src.pos(srcind(1:3),3),'r*');
axis equal

% for i=1:numel(L), %% double check indexing
%     tmpL=L{i};
%     Lx=tmpL(:,1:3:end);
%     Ly=tmpL(:,2:3:end);
%     Lz=tmpL(:,3:3:end);
%     L2=Lall{i};
%     rind=99;
%     L2{rind}-[Lx(:,rind) Ly(:,rind) Lz(:,rind)]
% end;
orstr=strvcat('lr','ap','si');
axstr={'X','Y','Z'}
Lind=[numel(L)];
maxall=zeros(length(Lind),3,3);
for m1=1:length(Lind),
    tmpL=L{Lind(m1)};
    %     Lx=tmpL(:,1:3:end);  % unpack order of 3 source orientations
    %     Ly=tmpL(:,2:3:end);
    %     Lz=tmpL(:,3:3:end);
    %
    maxLori=0;
    for dum=1:3,
        Lori=tmpL(:,dum:3:end); %% x, y or z orientation of source
        if max(max(abs(Lori)))>maxLori,
            maxLori=max(max(abs(Lori)));
        end;
    end;

    for s_ori=1:3, %% source orientation. x left right, y anterior posterior, z- superior inferior
        Lori=tmpL(:,s_ori:3:end); %% x, y or z orientation of source

        for ax1 = 1:3,

            ax=axstr(ax1);
            disp(ax{:})


            id = find(contains(grad.label,['-' ax{:}]));

            if ~isempty(id)

                comp2view = Lori(id,srcind);

                cfg = [];
                cfg.layout = lay;

                data.dimord    = 'chan_comp';
                data.topo      = comp2view;
                data.topolabel = {grad.label{id}}';
                data.time      = {1};

                cfg.component   = 1;
                cfg.interactive = 'no';
                cfg.comment     = 'no';
                cfg.title       = [models{Lind(m1)} sprintf('-%s-',orstr(s_ori,:)) ax{:}];

                % tmp_fig = figure('visible','off');
                figure
                maxall(m1,s_ori,ax1)=max(max(abs(data.topo)))
                if (maxall(m1,s_ori,ax1))==0, %% to stop plot crashing
                    warning('fixing data...')
                    models{m1}
                    data.topo(1,1)=1;
                end;
                % cfg.zlim=[-maxLori maxLori]/4;
                ft_topoplotIC(cfg,data);
                set(gcf,'color','w')
                %  colormap(brewermap(100,'RdBu'))
                colorbar
            end % isempty

        end % for ax
    end; s_ori
end; % for m1
