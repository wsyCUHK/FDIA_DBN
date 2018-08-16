function mincutatt
clc;
clear;

addpath('E:\Dropbox\MatlabCode\Graphics\TreePruning\matlab_bgl');
addpath('E:\Dropbox\MatlabCode\Graphics\MILP');

load('case14');

[m,n] = size(H);
[~,t]   = size(B);
Vs = n+1; %% the virtual sink

edgecost = [ones(1,9)*2,ones(1,11)]; 
edgecost(1)=10;
edgecost(2)=10;
AttackSet = [10,11];

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

cutEdge



