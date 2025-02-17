function src2 = tt_make_disk_sources(src, rad, nphi, phase)

% Make a ring of sources around the main cord
if nargin < 4
    phase = 0;
end

ang = linspace(0,(nphi-1)/nphi*2*pi,nphi) + phase;
x = rad * sin(ang);
y = rad * cos(ang);

disc = [x' y' zeros(size(y'))];

% Rather than let the disc be flat, we'll let it curve along the spine,
% like in tt_wrap_cord_v2, first plane of sources is flat
src2 = disc + src(1,:);
for ii = 2:length(src)
    vec = src(ii,:) - src(ii-1,:);
    vec = vec/norm(vec);
    orig_vec = [0 0 1];
    R = rodrigues(orig_vec,vec);
    rd = R*disc';
    pos = rd' + src(ii,:);
    src2 = [src2; pos];
end

end

function R = rodrigues(A, B)

A = A(:);
B = B(:);
v = cross(A, B);
c = dot(A, B);
vx = [  0    -v(3)  v(2);
    v(3)  0    -v(1);
    -v(2)  v(1)  0   ];
R = eye(3) + vx + vx^2 * (1 - c) / (norm(v)^2);

end

