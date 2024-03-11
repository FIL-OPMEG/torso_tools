function tt_add_bem
% add Helsinki BEM library to path;

addpath(fullfile(tt_path,'hbf_lc_p'));

if isempty(which('hbf_SetPaths'))
    disp('Cloning BEM library to repository!')
    !git submodule update --init
end

hbf_SetPaths;

% for FT compatibility install subfunctions to private folder for now
fnames = {'hbf_LFM_B_LC_xyz', 'hbf_Phiinf_xyz', 'hbf_Binf_xyz'};
exts = {'m','p'};
dir_in = fullfile(tt_path,'hbf_lc_p','hbf_calc','private');
dir_out = fullfile(tt_path,'private');
for ii = 1:numel(fnames)
    for jj = 1:numel(exts)
        fin = spm_file(fnames{ii},'ext',exts{jj},'path',dir_in);
        fout = spm_file(fin,'path',dir_out);
        copyfile(fin,fout);
    end
end