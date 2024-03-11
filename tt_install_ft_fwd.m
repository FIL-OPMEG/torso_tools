function tt_install_ft_fwd
% Adds forward solutions from FieldTrip to private folder if not there

funcs = sort({'fixpos','current_dipole','meg_leadfield1',...
    'surface_normals','fitsphere','meg_ini','legs',...
    'solid_angle','leadsphere_all','meg_forward'});
ft_private = fullfile(spm('dir'),'external',...
    'fieldtrip','forward','private');
tt_private = fullfile(tt_path,'private');

for ii = 1:numel(funcs)
    fprintf('Installing %s: ',funcs{ii});
    d = spm_select('fplist',ft_private,funcs{ii});
    if size(d,1)
        for jj = 1:size(d,1)
            try
                copyfile(deblank(d(jj,:)),tt_private);
            end
        end
        fprintf('OK\n');
    else
        fprintf('NOT FOUND\n');
    end
end

