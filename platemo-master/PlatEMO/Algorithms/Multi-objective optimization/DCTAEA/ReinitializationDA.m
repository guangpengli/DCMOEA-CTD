function Population = ReinitializationDA(Problem,Population)
% Re-initialize solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %至少有一个可行解,%50随机初始化
    
    evalues=0;
    
    R=Problem.Initialization(Problem.N/2);
    evalues=evalues+length(R);
    %选择以前的解50%
    Population=Problem.Evaluation(Population.decs);
    evalues=evalues+length(Population);
    
    [Population,~,~]   = EnvironmentalSelection([Population,R],Problem.N,2);
    Problem.FE=Problem.FE-evalues;
end