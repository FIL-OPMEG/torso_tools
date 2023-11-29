function [M]=tt_wrap_cord(sourcepos,rad)
%function wrap_cord(sourcepos,rad)
%% sourcespos are positions of sources (1 at each depth) moving along cord
%% rad is the radius of wrapper
%% M is the wrapper surface mesh

[dum,vertind]=max(std(sourcepos)); %% get main vertical axis

[dum,sortind]=sort(sourcepos(:,vertind));


cp=sourcepos(sortind,:); % in ascending order


astep=[0:pi/3:2*pi]
%astep=0;
allp=[];
for j=1:length(cp)-1,
    dsp=cp(j+1,:)-cp(j,:);
    magdsp=sqrt(dot(dsp,dsp));
    dsp=dsp./magdsp; %% unit vector pointing up cord
    dtan=cross(dsp,[0 1 0]); %% orthogonal to dsp
    dtan=dtan./sqrt(dot(dtan,dtan)); % unit
    drad=cross(dtan,dsp);
    for k=1:length(astep),
    
    p1=rad*drad*cos(astep(k))+rad*dtan*sin(astep(k))+cp(j+1,:)+dsp*magdsp*k/length(astep);
    allp=[allp; p1];

    end;
    
end;

k = boundary(allp(:,1),allp(:,2),allp(:,3),1);
M.faces=k;
M.vertices=allp;