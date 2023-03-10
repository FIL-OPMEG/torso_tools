function src = tt_generate_spine_grid(S)

if ~isfield(S,'subject'); error('please provide subject mesh!'); end
if ~isfield(S,'T'); error('please provide the transformation matrix!'); end
if ~isfield(S,'width');         S.width = 80;           end
if ~isfield(S,'depth');         S.depth = 10:10:80;     end
if ~isfield(S,'resolution');    S.resolution = 10;      end
if ~isfield(S,'mask');          S.mask = 1;             end

%-rotate the body scan for easier grid generation later
%-----------------------------------------------------
fids = tt_get_template_fids(S.T);

% halfway point between shoulder fids
hp = 0.5*(fids(2,:) + fids(1,:));
% vector between shoulder fids
hpv = fids(2,:) - fids(1,:);

% find the rotation in z which aligns the shoulders in y (is this a special
% case with just this scan?)
angD = -atand(hpv(2)/hpv(1));
R1 = rotmatZ(angD,hp);

% rotate the mesh for the grid generation
sub2 = spm_mesh_transform(S.subject,R1);
fids = tt_get_template_fids(R1*S.T);

% work out if the spine is in the +/-ve direction along y.
tmp = fids(3,:) - hp;
ydir = heaviside(tmp(2));

% get units to scale space with later;
unit = tt_determine_mesh_units(tt_load_meshes(R1*S.T));

%-Generate a plane of the spine, and raycast to get shape
%--------------------------------------------------------------------
bk = sub2.vertices;
x_midline = hp(1);
z_midline = 0.5*(min(bk(:,3))+max(bk(:,3)));
z_range = 0.8*range(bk(:,3));

% Generate plane
min_x = min([(x_midline - 0.5*S.width) (x_midline + 0.5*S.width)]);
max_x = max([(x_midline - 0.5*S.width) (x_midline + 0.5*S.width)]);
min_z = min([(z_midline - 0.5*z_range) (z_midline + 0.5*z_range)]);
max_z = max([(z_midline - 0.5*z_range) (z_midline + 0.5*z_range)]);

y_start = hp(2) + 10*tmp(2);
[xgrid,ygrid,zgrid] = meshgrid(min_x:10:max_x,y_start,min_z:10:max_z);
switch ydir
    case 1
        ray = [0 -1 0];
    case 0
        ray = [0 1 0];
end

G = gifti(sub2);
plane_x = [];
plane_y = [];
plane_z = [];
for ii = 1:numel(xgrid)
    R = struct('orig',[xgrid(ii) ygrid(ii) zgrid(ii)]',...
        'vec',ray');
    [I,P] = spm_mesh_ray_intersect(G,R);
    if ~isempty(P)
        if size(P,1) > 1
            switch ydir
                case 1
                    [~,pid] = max(P(:,2));
                case 0
                    [~,pid] = min(P(:,2));
            end
        else
            pid = 1;
        end
        plane_x(end+1) = P(pid,1);
        plane_y(end+1) = P(pid,2);
        plane_z(end+1) = P(pid,3);
    end
end

%-shift across depths
%--------------------
grid_all = [];

for ii = 1:numel(S.depth)
    tmp_plane = [plane_x; plane_y; plane_z]'+ ray*S.depth(ii);
    
    grid_all = cat(1,grid_all,tmp_plane);
    
end

%-mask sources to only be inside torso model
%-------------------------------------------
if S.mask
    fprintf('Determining sources inside torso model: ');
    [meshes,names] = tt_load_meshes(R1*S.T);
    mid = find(contains(names,'torso'));
    for ii = 1:length(grid_all)
        pos = grid_all(ii,:);
        inside(ii) = tt_is_inside(pos,...
            meshes{mid}.vertices,meshes{mid}.faces);
    end
    grid_all(~inside,:) = [];
    fprintf('COMPLETE\n');
end

%-make FT style source space container, transform back to orig space
%-------------------------------------------------------------------
src = [];
src.pos = grid_all;
src.inside = ones(length(grid_all),1);
src.unit = unit;
src = ft_transform_geometry(inv(R1),src);

end

function M = rotmatZ(deg,cp)
% make a rotation matrix around the Z axis centered at a point of choice
MT = [1 0 0 -cp(1);
    0 1 0 -cp(2);
    0 0 1 -cp(3);
    0 0 0 1];

MR = [cosd(deg) -sind(deg) 0 0;
    sind(deg) cosd(deg) 0 0
    0 0 1 0;
    0 0 0 1];

iMT = [1 0 0 cp(1);
    0 1 0 cp(2);
    0 0 1 cp(3);
    0 0 0 1];

M = iMT*MR*MT;
end

% bk=mesh.pos; %% back
% extrawidth=10;
% spxcentre=870; spwidth=25+extrawidth;
% spylim=-550;
% spzlim=250;
%
% plot3(bk(:,1),bk(:,2),bk(:,3),'go');
% xlabel('x');
% ylabel('y');
% zlabel('z');
%
% useind=intersect(find(bk(:,1)>(spxcentre-spwidth)),find(bk(:,1)<(spxcentre+spwidth)));
% useind=intersect(useind,find(bk(:,2)>spylim));
% useind=intersect(useind,find(bk(:,3)<spzlim));
%
% hold on
% plot3(bk(useind,1),bk(useind,2),bk(useind,3),'r.');
% halfplane=[bk(useind,1),bk(useind,2),bk(useind,3)];
% sf = fit([bk(useind,1), bk(useind,3)],bk(useind,2),'poly23');
% % quiver3(chanpos(:,1), chanpos(:,2),chanpos(:,3),...
% %     chanori(:,1), chanori(:,2), chanori(:,3),'color','b','linewidth',1)
%
%
% depthvals=-[10:10:70]; %% depth of plane to image onto
% pos = [];
%
% for depthind = 1:numel(depthvals)
%
%     coordoffset=depthvals(depthind); %% mm from back
%
%     gridstep=10;
%     lxind=0;lzind=0;
%     lzrange=min(bk(useind,3)):gridstep:max(bk(useind,3));
%     Nz=length(lzrange);
%     lxrange=min(bk(useind,1)):gridstep:max(bk(useind,1));
%     Nx=length(lxrange);
%     sourceind=zeros(Nx*Nz,2);
%     sourcepos=zeros(Nx*Nz,3);
%     count=0;
%     % figure;
%     pstep=1:10:length(bk);
%
%
%
%     for lzind=1:length(lzrange),
%         for lxind=1:length(lxrange)
%             count=count+1;
%             lx=lxrange(lxind);lz=lzrange(lzind);
%             %         plot3(lx,sf(lx,lz)+coordoffset,lz,'m.');
%             sourcepos(count,:)=[lx,sf(lx,lz)+coordoffset,lz];
%         end
%
%
%     end
%
%     pos = cat(1,pos,sourcepos);
%
% end
%
% plot3(pos(:,1),pos(:,2),pos(:,3),'k*');
%
% cfg = [];
% cfg.method = 'basedonpos';
% cfg.sourcemodel.pos = pos;
% cfg.unit = 'mm';
%
% src = ft_prepare_sourcemodel(cfg);
% src=ft_convert_units(src,'m');
%
% save('D:\thorax_model\sven\sourcespace_raw.mat');
