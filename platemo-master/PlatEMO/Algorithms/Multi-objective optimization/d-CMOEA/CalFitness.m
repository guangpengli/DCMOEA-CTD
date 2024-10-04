function [FrontNo,CrowdDis] = CalFitness(Population,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    zmin_o=min(Population.objs,[],1);
    zmax_o=max(Population.objs,[],1);
    PopObjs=(Population.objs-repmat(zmin_o,length(Population),1))./(repmat(zmax_o,length(Population),1)-repmat(zmin_o,length(Population),1));
    
    %约束初始化
    cons=max(0,Population.cons);
    maxCons=max(cons,[],1);
    cons=cons./(maxCons+0.000001);
    
    CV=sum(cons,2)/size(cons,2);
    zmin=min(CV);
    zmax=max(CV);
    CVs=(CV-repmat(zmin,length(Population),1))./(repmat(zmax,length(Population),1)-repmat(zmin,length(Population),1)+0.00001);
    fs=find(CV<=0);
    rf=length(fs)/length(Population);
    
    d=sqrt(CVs.^2+PopObjs.^2);
    
    Y=PopObjs;
    Y(fs,:)=0;
    
    p=(1-rf).*CVs+rf.*Y;
    F=d+p;

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(F,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(F,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end