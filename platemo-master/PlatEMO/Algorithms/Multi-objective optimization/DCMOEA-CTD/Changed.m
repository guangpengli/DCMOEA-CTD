function [changed,type] = Changed(Problem,Population)
% Detect whether the problem changes

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    RePop1  = Population(randperm(end,ceil(end/10)));
    
    RePop2Objs=Problem.CalObj(RePop1.decs);
    RePop2Cons=Problem.CalCon(RePop1.decs,RePop1.objs);
    

    
    changed=false;
    type=0;
    VO=~isequal(RePop1.objs,RePop2Objs);
    VC=~isequal(RePop1.cons,RePop2Cons);

    if  VO && ~VC
        % 1--基于目标分类

        
        RePopulation1Objs  = Problem.CalObj(Population.decs);
        [FrontNo,~] = NDSort(Population.objs,length(Population));  
        [FrontNo1,~] = NDSort(RePopulation1Objs,length(Population));
        
        if isequal(FrontNo,FrontNo1)
            %约束不变，UPS不变
            type=1;
        else
            %约束不变，UPS变
            type=4;
        end
        
    elseif ~VO && VC
        % 2
        type=2;
    elseif VO && VC
        % 3
        RePopulation1Objs  = Problem.CalObj(Population.decs);
        [FrontNo,~] = NDSort(Population.objs,length(Population));
        Fs=find(FrontNo==1);
        FrontNo=FrontNo(Fs);
        [FrontNo1,~] = NDSort(RePopulation1Objs(Fs,:),length(FrontNo));
         if isequal(FrontNo,FrontNo1)
            %约束变，UPS不变
            type=3;
        else
            %约束变，UPS变
            type=5;
        end
    end
    if type~=0
        changed=true;
    end
    
end