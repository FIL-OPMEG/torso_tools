function M = tt_register_torso(S)
% Generates the affine transform required to fit the canonical torso to a
% subject scan

if ~isfield(S,'fiducials'); error('you must provide fiducial locations'); end
if ~isfield(S,'subject'); error('you must provide a mesh to fit to!'); end
if ~isfield(S,'dist'); S.dist = 0.06; end
if ~isfield(S,'plot'); S.plot = false; end


% Load the canonical torso, generate fiducials
%----------------------------------------------
torso = export(gifti(fullfile(tt_path,'torso.gii')),'patch'); % Units: m

% fiducuals of the torso
% - left shoulder   (pt 901)
% - right shoulder  (pt 1953)
% - low spine       (pt 1147)
torso_fids = [torso.vertices(901,:);
    torso.vertices(1953,:);
    torso.vertices(1147,:)];

% Register part 0: get the two meshes into similar units
%--------------------------------------------------------
sf = determine_body_scan_units(S.fiducials,torso_fids);
M0 = [sf 0 0 0
    0 sf 0 0
    0 0 sf 0
    0 0 0 1];

torso = spm_mesh_transform(torso,M0);
% Update fiducials post scaling
torso_fids = [torso.vertices(901,:);
    torso.vertices(1953,:);
    torso.vertices(1147,:)];

if S.plot
    figure;clf
    patch(S.subject,'facecolor','none','edgecolor','k',...
        'edgealpha',0.3,'clipping','off');
    axis equal
    axis off
    set(gcf,'color','w');
    hold on
    h_fid   = plot3(S.fiducials(:,1),S.fiducials(:,2),...
        S.fiducials(:,3),'oc');
    set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');
end

% Register part 1: rigid body transform between fiducials
%--------------------------------------------------------
M1 = spm_eeg_inv_rigidreg(S.fiducials',torso_fids');
torso = spm_mesh_transform(torso,M1);
% Update fiducials post RBT
torso_fids = [torso.vertices(901,:);
    torso.vertices(1953,:);
    torso.vertices(1147,:)];

% Register part 2: ICP
%--------------------------------------------------------
[~,D] = knnsearch(S.subject.vertices,torso.vertices);
id = find(D<=S.dist*sf);
torso_fids = [torso.vertices(901,:);
    torso.vertices(1953,:);
    torso.vertices(1147,:)];

M2 = spm_eeg_inv_icp(S.subject.vertices',torso.vertices(id,:)',...
    S.fiducials',torso_fids',[],[],1);

if S.plot
    torso = spm_mesh_transform(torso, M2);
    torso_fids = [torso.vertices(901,:);
    torso.vertices(1953,:);
    torso.vertices(1147,:)];
    patch(torso,'facecolor','none','edgecolor','b','edgealpha',0.3)
    h_fidmr = plot3(torso_fids(:,1),torso_fids(:,2),...
        torso_fids(:,3),'dm');
    set(h_fidmr,'MarkerFaceColor','m',...
        'MarkerSize',12,'MarkerEdgeColor','k');
end

M = M2*M1*M0;

end

function sf = determine_body_scan_units(body_fids,torso_fids)
% Determine if the units of the torso and body are the same by looking at
% the triangle which is made between the fiducials

body_vec = (body_fids([1 2],:) - body_fids(3,:));
torso_vec = torso_fids([1 2],:) - torso_fids(3,:);
    
body_area = norm(cross(body_vec(1,:),body_vec(2,:)));
torso_area = norm(cross(torso_vec(1,:),torso_vec(2,:)));

pow = round(log10(sqrt(body_area/torso_area)));

sf = 10^pow;

end

