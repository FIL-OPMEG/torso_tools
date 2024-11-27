function R1 = tt_align_shoulders(T)

if ~nargin; error('please provide the transformation matrix!'); end

%-rotate the body scan for easier grid generation later
%-----------------------------------------------------
fids = tt_get_template_fids(T);

% halfway point between shoulder fids
hp = 0.5*(fids(2,:) + fids(1,:));
% vector between shoulder fids
hpv = fids(2,:) - fids(1,:);

% find the rotation in z which aligns the shoulders in y (is this a special
% case with just this scan?)
angD = -atand(hpv(2)/hpv(1));
R1 = rotmatZ(angD,hp);

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