clearvars
close all
clc

mesh = ft_read_headshape('D:\thorax_model\sven\seated_body_registered.stl');
% mesh = ft_convert_units(mesh,'m');

bk=mesh.pos; %% back
extrawidth=10;
spxcentre=870; spwidth=25+extrawidth;
spylim=-550;
spzlim=250;

plot3(bk(:,1),bk(:,2),bk(:,3),'go');
xlabel('x');
ylabel('y');
zlabel('z');

useind=intersect(find(bk(:,1)>(spxcentre-spwidth)),find(bk(:,1)<(spxcentre+spwidth)));
useind=intersect(useind,find(bk(:,2)>spylim));
useind=intersect(useind,find(bk(:,3)<spzlim));

hold on
plot3(bk(useind,1),bk(useind,2),bk(useind,3),'r.');
halfplane=[bk(useind,1),bk(useind,2),bk(useind,3)];
sf = fit([bk(useind,1), bk(useind,3)],bk(useind,2),'poly23');
% quiver3(chanpos(:,1), chanpos(:,2),chanpos(:,3),...
%     chanori(:,1), chanori(:,2), chanori(:,3),'color','b','linewidth',1)



depthvals=-[10:10:70]; %% depth of plane to image onto
pos = [];

for depthind = 1:numel(depthvals)
    
    coordoffset=depthvals(depthind); %% mm from back
    
    gridstep=10;
    lxind=0;lzind=0;
    lzrange=min(bk(useind,3)):gridstep:max(bk(useind,3));;
    Nz=length(lzrange);
    lxrange=min(bk(useind,1)):gridstep:max(bk(useind,1));
    Nx=length(lxrange);
    sourceind=zeros(Nx*Nz,2);
    sourcepos=zeros(Nx*Nz,3);
    count=0;
    % figure;
    pstep=1:10:length(bk);
    
    
    
    for lzind=1:length(lzrange),
        for lxind=1:length(lxrange)
            count=count+1;
            lx=lxrange(lxind);lz=lzrange(lzind);
            %         plot3(lx,sf(lx,lz)+coordoffset,lz,'m.');
            sourcepos(count,:)=[lx,sf(lx,lz)+coordoffset,lz];
        end
        
        
    end
    
    pos = cat(1,pos,sourcepos);
    
end

plot3(pos(:,1),pos(:,2),pos(:,3),'k*');

cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = pos;
cfg.unit = 'mm';

src = ft_prepare_sourcemodel(cfg);
src=ft_convert_units(src,'m');

save('D:\thorax_model\sven\sourcespace_raw.mat');
