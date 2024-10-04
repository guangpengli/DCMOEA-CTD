function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    cons=max(0,Population.cons);
    CV=sum(cons,2);
    zmin=min(CV,[],1);
    zmax=max(CV,[],1);
    CVs=(CV-repmat(zmin,length(Population),1))./(repmat(zmax,length(Population),1)-repmat(zmin,length(Population),1)+0.000001);
    
    fs=find(CV<=0);
    rf=length(fs)/N;
    zmin=min(Population.objs,[],1);
    zmax=max(Population.objs,[],1);
    PopObjs=(Population.objs-repmat(zmin,length(Population),1))./(repmat(zmax,length(Population),1)-repmat(zmin,length(Population),1));
    
    if rf==0
        d=CVs;
        X=0;
    else
        d=sqrt(CVs.^2+PopObjs.^2);
        X=CVs;
    end

    Y=PopObjs;
    Y(fs,:)=0;
    p=(1-rf).*X+rf.*Y;
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
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end