function GetParamGraph_RecG(Graph)

%getinfo about graph
nodes=Graph.nd;
edges=Graph.ed;

Holder=Graph.Data;
Fixed = Holder{1};
savdir = 'C:\Users\samue\Documents\MATLAB\networkgen\dataFiles';

% vector=Graph.rv;
N = length(nodes);
% Rx=vector(1);
% Ry=vector(2);
% Rz=vector(3);
coord1=nodes(:,1)';
coord2=nodes(:,2)';
coord3=nodes(:,3)';
T2=edges;
Nr=size(edges,1);

%% first part of table - coordinates of nodes
fileID = fopen('Settings.txt','r');
xT1=textscan(fileID,'%s','whitespace','\n');
fclose(fileID);
xT1=xT1{1,1};

xT1=[cell(1,0)
'<?xml version="1.0" encoding="UTF-8"?>'
'<data>'
xT1
'    <nodes default-mass="1.0">'];

% xT1=textscan('Setings.txt','%s');
for i=1:N
    xT1=[xT1;['        <node>',num2str(coord1(i)),' ',num2str(coord2(i)),' ',num2str(coord3(i)),'</node>']];
end

%%xml file
T1=[cell(1,0),'nodes'];

% if max(max(nodes))>1000
    
% else
%     coeff=1;
% end

for i=1:N
    T1=[T1;[num2str(i),' x =   ',num2str(coord1(i)),'     y =    ',num2str(coord2(i)),'     z =    ',num2str(coord3(i))]];
end

%% second part of table - numbers of nodes for all links
T1=[T1;cell(1,1);'elements of type Beam_Spring_2 '];
for i=1:size(T2,1)
    T1=[T1;[num2str(i),' nodes =  [    ',num2str(T2(i,1)),',    ',num2str(T2(i,2)),'] material = Material1']];
end

%%for xml file
xT1=[xT1
'    </nodes>'
'    <links>'];
for i=1:size(T2,1)
    xT1=[xT1;['        <link>',num2str(T2(i,1)-1),' ',num2str(T2(i,2)-1),'</link>']];
end
xT1=[xT1;['    </links>']];
% [~,mx]=max([Rx,Ry,Rz]);mx=mx(1);
% 10 percent threshold

%notice we sort on the x axis so we have force vector [-1,0,0]
all=[coord1',coord2',coord3'];
%in2 = find(all(:,1)<=max(all(:,1))*Fixed);

intForce = [];
targets = [];



%For the fixed points we choose the ones on the bottom of the square
xT1=[xT1  
'    <fixed>'];
%for i=1:length(in2)
 %   xT1=[xT1;['            <node>',num2str(in2(i)-1),'</node>']];
%end

xT1=[xT1
'    </fixed>'
'</data>'];

lay=cell(1,0);
for i=1:size(T2,1)
    lay=[lay;[num2str(T2(i,1)-1),' ',num2str(T2(i,2)-1)]];
end
coeff=1000/max(max(nodes));
for i=1:size(nodes,1)
    lay=[lay;['//NODECOORD ',num2str(i-1),' ',num2str(round(coeff*nodes(i,1))),' ',num2str(round(coeff*nodes(i,2))),' ',num2str(round(coeff*nodes(i,3)))]];
end
%% write to files
fileID = fopen(fullfile(savdir,['Graph_',Graph.String,'.txt']),'w');
C=T1;
[nrows,~] = size(C);
formatSpec = '%s\n';
for row = 1:nrows
    fprintf(fileID,formatSpec,C{row,:});
end

fileID = fopen(fullfile(savdir,['data_',Graph.String,'_4.xml']),'w');
C=xT1;
[nrows,~] = size(C);
formatSpec = '%s\n';
for row = 1:nrows
    fprintf(fileID,formatSpec,C{row,:});
end

% fileID = fopen(fullfile(savdir,['BioGraph_',Graph.String,'.layout']),'w');
% C=lay;
% [nrows,~] = size(C);
% formatSpec = '%s\n';
% for row = 1:nrows
%     fprintf(fileID,formatSpec,C{row,:});
% end


%% compute degree distribution
col=[];
nod=[];
if sum(edges)
    nod=[edges(:,1);edges(:,2)];
end
for i=1:N
    col=[col;length(find(nod==i))];
end
un=unique(col);
if min(un)>0
    un=[0;un];
end
distr=zeros(1,max(un));
for i=un'
    distr(i+1)=length(find(col==i))/N;
end

%% plot degree distribution
figure


if isfield(Graph,'vr')
    hold on
    varr=Graph.vr;
    ro=varr(1);
    a=varr(2);
    expmu=varr(3);
    sigma=varr(4);
    NumIntLen=varr(5);
    MaxDeg=varr(6);
    Rad=varr(7);
    Len=varr(8);
    VectEd=3:MaxDeg;
    RaspEd=VectEd.^-a;
%     RaspEd([1,2,3])=0;
    if RaspEd(1)==0
        RaspEd(1)=1;
    end
    RaspEd=RaspEd/sum(RaspEd);
    plot(VectEd,RaspEd,'r','LineWidth',2)
end

bar(0:max(un),distr);
title('Nodes degree distribution');
xlabel('Number of Neighbors');
ylabel('Normalized Number of Nodes');


if isfield(Graph,'vr')
    plot(VectEd,RaspEd,'r','LineWidth',2)
    legend(['Power law \alpha =',num2str(a)],'Generated network');
end

%% compute length distribution
for i=1:size(T2,1)
    len(i)=sqrt((coord1(T2(i,1))-coord1(T2(i,2)))^2+(coord2(T2(i,1))-coord2(T2(i,2)))^2+(coord3(T2(i,1))-coord3(T2(i,2)))^2);
end

%% plot length distribution
figure
% [p,w]=hist(len,ceil(max(len)));


if isfield(Graph,'vr')
    hold on
    xmax=sqrt(Len^2+4*Rad^2);
% 	IntLen=xmax/NumIntLen;
    xx=0:xmax/5000:xmax;
    y=lognpdf(xx,log(expmu),sigma);
%     xx=IntLen/2:IntLen:xmax-IntLen/2;
%     yy=lognpdf(xx,log(expmu),sigma);
%     y=y/sum(yy);
    plot(xx,y,'r','LineWidth',2);

end




if isfield(Graph,'ni')
    NumIntLen=Graph.ni;
    xmax=sqrt(Len^2+4*Rad^2);
    IntLen=xmax/NumIntLen;
    x=IntLen/2:IntLen:xmax-IntLen/2;
    [p,w]=hist(len,x);
else
    [p,w]=hist(len,11);
end
bar(w,p/size(T2,1)/(w(2)-w(1)));

set(gca, 'FontSize', 40)
title('Edges length distribution');
xlabel('Length');
ylabel('Number of Edges');


if isfield(Graph,'vr')
    plot(xx,y,'r','LineWidth',2);
    legend(['Lognormal \mu=ln(',num2str(expmu),') \sigma=',num2str(sigma)],'Generated network');
end
%% compute volume
% [~,mx]=max(vector);mx=mx(1);
% xyz=[1,2,3];
% mn=xyz(xyz~=mx);
% vector=max(nodes);
% central=vector(mn)/2;
% for  i=1:N
%     tt=nodes(i,:);
%     ttmn=tt(mn);
%     rads(i)=sqrt((ttmn(1)-central(1))^2+(ttmn(2)-central(2))^2);
% end
DateString=Graph.String;



