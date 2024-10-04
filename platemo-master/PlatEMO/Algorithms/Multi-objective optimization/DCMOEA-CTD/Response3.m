function [Population1,Population2,evalues,Fitness1,Fitness2]=Response3(Pop1,Pop2,Problem,UPS,t)
    %% Type4 约束函数固定，UPS变化，Pop2差分预测，Pop1多样性策略
    %Pop1
    zmin1=min(Pop1.decs,[],1);
    zmax1=max(Pop1.decs,[],1);
    Gauss1=normrnd(0,1,[Problem.N/2,1]);
    index=randperm(Problem.N,50);
    X1=Pop1(index).decs+Gauss1.*(zmax1-zmin1);
    P1=Problem.Initialization(Problem.N/2);
    Population1=Problem.Evaluation(X1);

    
    %Pop2
    %策略1
    Pop2t_1=UPS(t-1,:);
    [FrontNo,~] = NDSort(Pop2t_1.objs,length(Pop2t_1));
    NS=Pop2t_1(FrontNo==1);
    C2t_1=mean(NS.decs,1);
    
    [FrontNo,~] = NDSort(Pop2.objs,length(Pop2));
    NS=Pop2(FrontNo==1);
    C2t=mean(NS.decs,1);

    Delta_t=C2t-C2t_1;
    X2=Pop2.decs;

    X2=X2+repmat(Delta_t,size(X2,1),1);
    Population2=Problem.Evaluation(X2);
    
    evalues=200;
    [Population1,Fitness1]=EnvironmentalSelection([Population1,Population2,P1],Problem.N,1);
    [Population2,Fitness2]=EnvironmentalSelection([Population1,Population2,P1],Problem.N,2);
   
   
 
end