classdef LPFEMMD<ALGORITHM
 % <multi> <real/binary/permutation> <constrained/none>   
 
 methods
     function main(Algorithm,Problem)
         %% Parameter setting
        K = Algorithm.ParameterSet(20);
         
        %% Generate random population and weight vectors
       
        [W,K]=UniformPoint(K,Problem.M);
        %ceil向上取整，S是每个子区域的个数
        Problem.N  = ceil(Problem.N/K)*K;
        S          = Problem.N/K;
        Population = Problem.Initialization();
        a=true;
        p=zeros(1,K);
        gen=0;
        %Fitness    = CalFitness(Population.objs);
        %% Optimization
        while Algorithm.NotTerminated(Population)
            gen=gen+1;
            MatingPool=randi(length(Population),1,Problem.N);
            Offspring=OperatorGA(Population(MatingPool));
%             MatingPool = TournamentSelection(2,Problem.N,Fitness);
%             Offspring  = OperatorGA(Population(MatingPool));
%             MatingPoolLocal      = randi(S,S,K) + repmat(0:S:S*(K-1),S,1);
%             MatingPoolGlobal     = randi(Problem.N,1,Problem.N);
%             rnd                  = rand(S,K) < 0.7;
%             MatingPoolLocal(rnd) = MatingPoolGlobal(rnd);
%             Offspring  = Operator(Population,Population(MatingPoolLocal(:)));
            Q=[Population,Offspring];
            
            %标准化目标空间
            zmin=min(Q.objs,[],1);
            zmax=max(Q.objs,[],1);
            qobjs=(Q.objs-repmat(zmin,length(Q),1))./(repmat(zmax,length(Q),1)-repmat(zmin,length(Q),1));
            flag=false;
            if mod(gen,20) == 0||a
                flag=true;
                a=false;
            end
            
            %将Q划分为多个区域
            [Population,p]=Associate(Q,qobjs,W,S,flag,p);
            
        end
     end
 end
end