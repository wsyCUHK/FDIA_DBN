clear;    
    
addpath('E:\Dropbox\MatlabCode\Graphics\TreePruning\matlab_bgl');
addpath('E:\Dropbox\MatlabCode\Graphics\MILP');
addpath('E:\Dropbox\MatlabCode\matpower');
labelAtt=zeros(10000,180);
StateRe=zeros(10000,180);

netsize=118;
for datai=1:100000
mpc = loadcase('case118.m');
    %[J, Ybus, Yf, Yt] = makeJac(mpc);
    
    
    
    
    for k=1:netsize
        rat=0.5*rand+1;
        %j=randi([2 14]);
        %realreact=randi([3 4]);
        mpc.bus(k, 3) = mpc.bus(k, 3)*rat;
    end
    results=rundcopf(mpc);
    %injectionpower(i,:)=results.branch(:,14);
    AA=full(makeBdc(mpc));
    [J, Ybus, Yf, Yt] = makeJac(mpc);
    B=imag(full(Ybus));
    

    

 
    load('E:\Dropbox\MatlabCode\Graphics\Topology\case118\case118.mat')
    
     CCC=zeros(118,118);
for i=1:186
for j=1:118
if A(i,j)>0.1
index1=j;
elseif A(i,j)<-0.1
index2=j;
end
end
CCC(index1,index2)=CCC(index1,index2)+1;CCC(index2,index1)=CCC(index2,index1)+1;
end
         wah=zeros(118,118);
wah(find(CCC>1))=1;
    
    DDmatrix=zeros(20,1);
    indexD=1;
    for conD=1:size(AA,2)
        for conDD=(conD+1):size(AA,1)
            if AA(conD,conDD) ~=0
                if wah(conD,conDD)==0
                DDmatrix(indexD)=AA(conD,conDD);
                indexD=indexD+1;
                else
                DDmatrix(indexD)=AA(conD,conDD);
                indexD=indexD+1;
                DDmatrix(indexD)=AA(conD,conDD);
                indexD=indexD+1;
                end
            end
        end
    end
    D=diag(DDmatrix);
    
    
    PB=D*A;
    
    MatrixSizeH=length(AA);
    for ii=1:MatrixSizeH
        for j=1:MatrixSizeH
            if ii==j
                H1(ii,j)=sum(AA(ii,:))-AA(ii,ii);
                %H2(ii,j)=sum(B(ii,:))-B(ii,ii);
            else
                H1(ii,j)=-AA(ii,j);
                %H2(ii,j)=-B(ii,j);
            end
        end
    end
    BigH=[H1;PB];
    
    %Hsize=length(B,1);
    %SelectedH=zeros(Hsize,14);
    load('E:\Dropbox\MatlabCode\Graphics\Topology\case118\case118.mat')
    FlS = FlowSet(B,InjetIdx); % The set of branches with measurements
    BS=[]; %The set of buses with measurements
    [m6 n6]=size(B);
   for m66=InjetIdx:m6
   TempIndex=find(B(m66,:));
   [doesmatter numofbus]=max(sum(abs(A(TempIndex,:))));
   BS(m66-InjetIdx+1)=numofbus;
   end
       
%     for SeH=1:Hsize
%     SelectedH(SeH,:)=BigH(,:)
%     end
SelectedRow=[BS FlS+14];
     SelectedH1=BigH(SelectedRow,:);
%     SelectedH2=BigH2(SelectedRow,:);
%     SelectedH3=BigH3(SelectedRow,:);
%     SelectedH4=BigH4(SelectedRow,:);
    
    

% RealStateTran=zeros(length(results.branch(:,1)),1);
%  RealStateTran(:)=results.branch(:,14);
    
    RealStateTran=zeros(length(results.branch(:,1)),3);
    RealState=zeros(118,1);
    RealStateTran(:,1)=results.branch(:,1);
    RealStateTran(:,2)=results.branch(:,2);
    RealStateTran(:,3)=results.branch(:,14);
   AAnglePhase(datai,:)=results.bus(:,9);
    
    for Reali=1:118
        RealState(Reali)=sum(RealStateTran(find(RealStateTran(:,1)==Reali),3))-sum(RealStateTran(find(RealStateTran(:,2)==Reali),3));
    end
    %injectionpower(i,:)=RealState;
    State=[RealState(BS);-RealStateTran(FlS,3)];
    
    
    
    v=ones(19,1);
    v=2*v;
    VV=diag(v);
  %  BBB1(1)=norm(State-SelectedH1*pinv(SelectedH1'*VV*SelectedH1)*SelectedH1'*VV*State);
%       BBB2(i)=norm(State-SelectedH2*pinv(SelectedH2'*VV*SelectedH2)*SelectedH2'*VV*State);
%         BBB3(i)=norm(State-SelectedH3*pinv(SelectedH3'*VV*SelectedH3)*SelectedH3'*VV*State);
%           BBB4(i)=norm(State-SelectedH3*pinv(SelectedH3'*VV*SelectedH3)*SelectedH3'*VV*State);
load('case118');

[m,n] = size(H);
[~,t]   = size(B);
Vs = n+1; %% the virtual sink



branchnum=size(A,1);
edgecost = [ones(1,floor(branchnum/2))*2,ones(1,branchnum-floor(branchnum/2))]; 
edgecost(1)=10;
edgecost(2)=10;
edgecost(floor(branchnum/4))=5;
edgecost(floor(branchnum/4)+1)=5;


AttackSet1 = [randi(netsize-2)+2,randi(netsize-2)+2,randi(netsize-2)+2,randi(netsize-2)+2,randi(netsize-2)+2,randi(netsize-2)+2,randi(netsize-2)+2];
AttackSet=unique(AttackSet1);

%% The set of lines with flow measurements
FlS = FlowSet(B,InjetIdx);

%% The set of edges measured by flow measurements
if m==InjetIdx
    InS = find(B(InjetIdx,:)~=0);
else
    InS=find(sum(B(InjetIdx:m,:))~=0);
end

MLine= setminus(1:t, setminus(1:t,union(FlS,InS)));%% indices of measured lines

%%the effective bus-to-bus incidence matrix C (TREE), removing redundant lines
C=RemainInd(MLine,A,edgecost);
C1 = extendGraph(C,AttackSet);

Z=sparse(C1);

[flowval cut R F] = max_flow(Z,1,Vs);
sourceSide =[];
sinkSide=[];

for i=1:n %% only the first n effective nodes
    if cut(i)==1
        sourceSide=union(i,sourceSide);
    else
        sinkSide=union(i,sinkSide);
    end
end

cutEdge=[];

for i=1:length(MLine)
    p = MLine(i);
    p1 = find(A(MLine(i),:)); %% find the ends of a line
    if (ismember(p1(1),sourceSide) && ismember(p1(2),sinkSide) ) || (ismember(p1(2),sourceSide) && ismember(p1(1),sinkSide) )
        cutEdge = union(p,cutEdge);
    end
end

  attackc=zeros(netsize,1);
  attackc(AttackSet)=1;
  attacka=SelectedH1*attackc;
  
  StateRe(datai,:)=State+attacka;
  labelAtt(datai,find(attacka))=1;
  datai
end
for datai=100001:110000
mpc = loadcase('case118.m');
    %[J, Ybus, Yf, Yt] = makeJac(mpc);
    
    
    
    
    for k=1:netsize
        rat=0.5*rand+1;
        %j=randi([2 14]);
        %realreact=randi([3 4]);
        mpc.bus(k, 3) = mpc.bus(k, 3)*rat;
    end
    results=rundcopf(mpc);
   
    
    RealStateTran=zeros(length(results.branch(:,1)),3);
    RealState=zeros(netsize,1);
    RealStateTran(:,1)=results.branch(:,1);
    RealStateTran(:,2)=results.branch(:,2);
    RealStateTran(:,3)=results.branch(:,14);
   AAnglePhase(datai,:)=results.bus(:,9);
    
    for Reali=1:netsize
        RealState(Reali)=sum(RealStateTran(find(RealStateTran(:,1)==Reali),3))-sum(RealStateTran(find(RealStateTran(:,2)==Reali),3));
    end
    %injectionpower(i,:)=RealState;
    State=[RealState(BS);-RealStateTran(FlS,3)];
    
  
  
  StateRe(datai,:)=State;
datai
end
labelAtt(100001:110000,:)=zeros(10000,180);

y2=wgn(110000,180,-10);
StateRe=StateRe+y2;
index_label=randperm(110000);
y_train=labelAtt(index_label(1:100000),:);
y_test=labelAtt(index_label(100001:110000),:);
x_train=StateRe(index_label(1:100000),:);
x_test=StateRe(index_label(100001:110000),:);
 
 filename = 'E:\Dropbox\Python\Data\data118_n0db.mat';
 save(filename,'x_train','x_test','y_test','y_train');
 
 y2=wgn(110000,180,-10);
StateRe=StateRe+y2;
index_label=randperm(110000);
y_train=labelAtt(index_label(1:100000),:);
y_test=labelAtt(index_label(100001:110000),:);
x_train=StateRe(index_label(1:100000),:);
x_test=StateRe(index_label(100001:110000),:);
 
 filename = 'E:\Dropbox\Python\Data\data118_n1db.mat';
 save(filename,'x_train','x_test','y_test','y_train');
 
  y2=wgn(110000,180,-10);
StateRe=StateRe+y2;
index_label=randperm(110000);
y_train=labelAtt(index_label(1:100000),:);
y_test=labelAtt(index_label(100001:110000),:);
x_train=StateRe(index_label(1:100000),:);
x_test=StateRe(index_label(100001:110000),:);
 
 filename = 'E:\Dropbox\Python\Data\data118_n2db.mat';
 save(filename,'x_train','x_test','y_test','y_train');

   y2=wgn(110000,180,-10);
StateRe=StateRe+y2;
index_label=randperm(110000);
y_train=labelAtt(index_label(1:100000),:);
y_test=labelAtt(index_label(100001:110000),:);
x_train=StateRe(index_label(1:100000),:);
x_test=StateRe(index_label(100001:110000),:);
 
 filename = 'E:\Dropbox\Python\Data\data118_n3db.mat';
 save(filename,'x_train','x_test','y_test','y_train');
    y2=wgn(110000,180,-10);
StateRe=StateRe+y2;
index_label=randperm(110000);
y_train=labelAtt(index_label(1:100000),:);
y_test=labelAtt(index_label(100001:110000),:);
x_train=StateRe(index_label(1:100000),:);
x_test=StateRe(index_label(100001:110000),:);
 
 filename = 'E:\Dropbox\Python\Data\data118_n4db.mat';
 save(filename,'x_train','x_test','y_test','y_train');