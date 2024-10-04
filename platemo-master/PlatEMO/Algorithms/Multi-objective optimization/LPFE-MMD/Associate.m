function [Population,p]=Associate(Population,pobjs,W,S,flag,p)
    K=size(W,1);
    %% Allocation of solutions to subproblems
    [~,transformation]=max(1-pdist2(pobjs,W,'cosine'),[],2);
    partition=zeros(S,K);
    for i=1:K
        current=find(transformation==i);
        % 第i个子域的解的个数不够S个
        if length(current)<S
            current = [current;randi(length(Population),S-length(current),1)];
        % 第i个子域的解个数超过S个
        elseif length(current)>S
            %计算当前子域的曲率信息，选择指标函数，选择下一代个体。
            %用当前子域中的非支配解去预测
            [FrontNo,~]=NDSort(pobjs(current),1);
            front1=find(FrontNo==1);
            
            front1objs=pobjs(current,:);
            front1objs(front1);
            
            if flag
                pi=CalculateCurvature(front1objs(front1,:));
                p(i)=pi;
            end
            
            if p(i)>1
                % 凹形
                con=max(pobjs(current),[],2);
            else
                % 凸形
                %收敛性
                con=sum(pobjs(current),2);
            end
            %分布性=1/（最近k邻的距离+2）
            distance=pdist2(pobjs(current),pobjs(current));
            distance(logical(eye(length(distance)))) = inf;
            distance = sort(distance,2);
            div = 1./(distance(:,floor(sqrt(length(current))))+2);
            fitness=con+div;
            [~,rank]=sort(fitness);
            current=current(rank(1:S));
        end
        partition(:,i)=current;
    end
    Population=Population(partition(:));
    
end