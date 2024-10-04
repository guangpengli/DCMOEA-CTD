function Population = Reinitialization(Problem,Population,Z)
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
    Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population(Next).objs,Population(Last).objs,Problem.N/2-sum(Next),Z,Zmin);
    Next(Last(Choose)) = true;
    % Population for next generation
    CS = Population(Next);
    NewPop=[NewPop CS];
    NewCVs=sum(max(0,NewPop.cons),2);
    newfs=find(NewCVs<=0);
    NewFs=NewPop(newfs);
    decs=NewPop.decs;
    zmin=min(decs,[],1);
    zmax=max(decs,[],1);
    Offsprings=[];
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
    NewPop=[NewPop Offsprings];
    
    Zmin          = min(NewPop(all(NewPop.cons<=0,2)).objs,[],1);
    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end
    Population   = EnvironmentalSelection(NewPop,length(Population),Z,Zmin);
    Problem.FE=Problem.FE-evalues;
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end