for datai=104001:110000
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
end
labelAtt(104001:110000,:)=zeros(6000,180);
y2=
StateRe=StateRe+y2;
index_label=randperm(110000);
y_train=labelAtt(index_label(1:100000),:);
y_test=labelAtt(index_label(100001:110000),:);
x_train=StateRe(index_label(1:100000),:);
x_test=StateRe(index_label(100001:110000),:);
 
 filename = 'E:\Dropbox\Python\Data\data118_ns.mat';
 save(filename,'x_train','x_test','y_test','y_train');