function tt_add_bem
% add Helsinki BEM library to path;

addpath(fullfile(tt_path,'hbf_lc_p'));

if isempty(which('hbf_SetPaths'))
    disp('Cloning BEM library to repository!')
    !git submodule update --init
end

hbf_SetPaths;