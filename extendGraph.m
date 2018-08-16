function C1=extendGraph(C,AttackSet)

infty=10000;
[n,~]=size(C);
C1 = zeros(n+1,n+1);

C1(1:n,1:n)=C;

for i=1:length(AttackSet)
    p = AttackSet(i);
    C1(p,n+1)=infty;
    C1(n+1,p)=infty;
end