function [M]=tt_wrap_cord_v2(sourcepos, rad, subdiv, Nphi)
% function wrap_cord(sourcepos,rad)
% sourcespos are positions of sources (1 at each depth) moving along cord
% rad is the radius of wrapper
% M is the wrapper surface mesh
if nargin < 3
    subdiv = 1;
    Nphi = 20;
elseif nargin < 4
    Nphi = 20;
end

% Create phantom halfway points between sources
if subdiv > 0
    for ss = 1:subdiv
        sourcepos2 = sourcepos(1,:);
        for ii = 2:length(sourcepos)
            half = mean(sourcepos(ii-1:ii,:));
            sourcepos2 = [sourcepos2;half;sourcepos(ii,:)];
        end
        sourcepos = sourcepos2;
    end
end

% sourcepos = sourcepos(1:10,:); % dont uncomment this - just for testing!

% Make an almost cylindrical tube
[z,phitemp]=meshgrid(sourcepos(:,3),linspace(0,(Nphi-1)/Nphi*2*pi,Nphi));
phitemp=phitemp(:);z=z(:);
etemp=delaunay(phitemp,z);

% Add the missing faces to close the tube
for ii = 1:length(sourcepos)-1
    tri_id = [Nphi*(ii-1)+1 Nphi*(ii+1), Nphi*ii;
        Nphi*(ii-1)+1 Nphi*ii+1 Nphi*(ii+1)];
    etemp = [etemp; tri_id];
end


xshift = repmat(sourcepos(:,1),1,Nphi)';
yshift = repmat(sourcepos(:,2),1,Nphi)';
x = xshift(:) + rad*cos(phitemp);
y = yshift(:) + rad*sin(phitemp);


% Currently this shears the tube, lets try and rotate each disk
xmat = reshape(x,Nphi,[]);
ymat = reshape(y,Nphi,[]);
zmat = reshape(z,Nphi,[]);

for ii = 2:length(sourcepos)
    pos = [xmat(:,ii) ymat(:,ii) zmat(:,ii)] - sourcepos(ii,:);
    vec = sourcepos(ii,:) - sourcepos(ii-1,:);
    vec = vec/norm(vec);
    orig_vec = [0 0 1];
    R = rodrigues(orig_vec,vec);
    pos = R*pos';
    pos = pos + sourcepos(ii,:)';
    xmat(:,ii) = pos(1,:);
    ymat(:,ii) = pos(2,:);
    zmat(:,ii) = pos(3,:);
end
x = xmat(:);
y = ymat(:);
z = zmat(:);

cyl.faces = etemp;
cyl.vertices = [x y z];


% cap each end with an isosphere?
hmesh = hemisphere_cap(rad, Nphi);
vid = find(abs(hmesh.vertices(:,3)) < 1e-3);
hmesh_ring = hmesh.vertices(vid,:);

% bottom cap first
% register the cap to the main spine
bottom_ring = flip([xmat(:,1) ymat(:,1) zmat(:,1)]);
M1 = spm_eeg_inv_rigidreg(bottom_ring',hmesh_ring');
bottom_cap = spm_mesh_transform(hmesh,M1);
% work out which points are replaced with spine points
hmesh_reg_ring = bottom_cap.vertices(vid,:);
lut = knnsearch(cyl.vertices,hmesh_reg_ring);
% append to lists
nverts = length(cyl.vertices);
bottom_cap.vertices(vid,:) = [];
bottom_cap.faces = bottom_cap.faces + nverts;
for ii = 1:numel(vid)
    id = find(bottom_cap.faces == vid(ii)+nverts);
    bottom_cap.faces(id) = lut(ii);
end
cyl.vertices = [cyl.vertices; bottom_cap.vertices];
cyl.faces = [cyl.faces; bottom_cap.faces];

% repeat for top!
top_ring = [xmat(:,end) ymat(:,end) zmat(:,end)];
M2 = spm_eeg_inv_rigidreg(top_ring',hmesh_ring');
top_cap = spm_mesh_transform(hmesh,M2);
hmesh_reg_ring = top_cap.vertices(vid,:);
lut = knnsearch(cyl.vertices,hmesh_reg_ring);
nverts = length(cyl.vertices);
top_cap.vertices(vid,:) = [];
top_cap.faces = top_cap.faces + nverts;
for ii = 1:numel(vid);
    id = find(top_cap.faces == vid(ii)+nverts);
    top_cap.faces(id) = lut(ii);
end
cyl.vertices = [cyl.vertices; top_cap.vertices];
cyl.faces = [cyl.faces; top_cap.faces];
%
% [M.vertices, M.faces] = hbf_CorrectTriangleOrientation(cyl.vertices, cyl.faces);
M = cyl;

end
function R = rodrigues(A, B)

if A == B

    R = eye(3);

else

    A = A(:);
    B = B(:);
    v = cross(A, B);
    c = dot(A, B);
    vx = [  0    -v(3)  v(2);
        v(3)  0    -v(1);
        -v(2)  v(1)  0   ];
    R = eye(3) + vx + vx^2 * (1 - c) / (norm(v)^2);
end
end

function hmesh = hemisphere_cap(rad, Nphi)

theta = linspace(0, (Nphi-1)/Nphi*2*pi, Nphi); % Azimuthal angle
phi = linspace(1/5*pi/2, pi/2, 5);   % Polar angle

% Create grid of points
[theta, phi] = meshgrid(theta, phi);
theta = theta';
theta = theta(:);
phi = phi';
phi = phi(:);
e = delaunay(theta,phi);
% Add the missing faces to close the walls
for ii = 1:4
    tri_id = [Nphi*(ii-1)+1 Nphi*(ii+1), Nphi*ii;
        Nphi*(ii-1)+1 Nphi*ii+1 Nphi*(ii+1)];
    e = [e; tri_id];
end
% rearrange so triangles are oriented correctly
e = e(:,[1 3 2]);

% Convert to Cartesian coordinates
x = rad * sin(phi) .* cos(theta);
y = rad * sin(phi) .* sin(theta);
z = rad * cos(phi);

% add the apex point
x = [0; x];
y = [0; y];
z = [rad; z];
% close the lid
e = e + 1;
l = 2:Nphi+1;
for ii = 1:Nphi
    tri_id = [1 l(1) l(2)];
    e = [e; tri_id];
    l = circshift(l,-1);
end

hmesh.vertices = [x y z];
hmesh.faces = e;

end