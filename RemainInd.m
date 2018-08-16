function C=RemainInd(MLine,A,edgecost)
%% construct a bus-bus incidence matrix D, bidirected
t = length(MLine);
[~,n]= size(A);
T = zeros(n,n); 

for i=1:t
    c = edgecost(MLine(i)); %% the cost of this line
    temp = find(A(MLine(i),:));    
    T(temp(1),temp(2))=c;
end

C=max(T,T');

