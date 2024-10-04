function Population = ReinitializationCA(Problem,Population)
% Re-initialize solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    flag=true;
    %至少有一个可行解,%50随机初始化
    NewPop=[];
    n=1;
    while flag
        R=Problem.Initialization(Problem.N/2);
        Problem.FE=Problem.FE-length(R);
        NewCVs=sum(max(0,R.cons),2);
        fs=find(NewCVs<=0);
        if ~isempty(fs)||n>10
            flag=false;
        end
        n=n+1;
    end
    evalues=0;
    NewPop=[NewPop R];
    p=length(fs)/(Problem.N/2);
    
    %选择以前的解50%
    Population=Problem.Evaluation(Population.decs);
    evalues=evalues+length(Population);
    PopObjs=Population.objs;
    [FrontNo,MaxFNo] = NDSort(PopObjs,Population.cons,Problem.N/2);
    Next = FrontNo < MaxFNo;
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(PopObjs,FrontNo);
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:Problem.N/2-sum(Next)))) = true;
    CS=Population(Next);
    NewPop=[NewPop CS];
    NewCVs=sum(max(0,NewPop.cons),2);
    newfs=find(NewCVs<=0);
    NewFs=NewPop(newfs);
    decs=NewPop.decs;
    zmin=min(decs,[],1);
    zmax=max(decs,[],1);
    Offsprings=[];
    if ~isempty(newfs)
        for i=1:length(CS)
            indx=randperm(length(newfs),1);
            xbest=NewFs(indx);
            xi=CS(i);
            xi_dec=xi.dec;
            xbest_dec=xbest.dec;
            for j=1:Problem.D
                xij=xi_dec(j)+rand*(xbest_dec(j)-xi_dec(j));
                r=rand*p;
                if xij<zmin(j)
                    xij=zmin(j)+r*(xij-zmin(j));
                elseif xij>zmax(j)
                    xij=zmax(j)+r*(xij-zmax(j));
                end
                xi_dec(j)=xij;
            end
            offspring=Problem.Evaluation(xi_dec);
            Offsprings=[Offsprings offspring];
        end
    end
    evalues=evalues+length(Offsprings);
    [Population,~,~]   = EnvironmentalSelection([NewPop,Offsprings],length(Population),1);
    Problem.FE=Problem.FE-evalues;
end