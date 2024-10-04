classdef DCTAEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
% Two-archive evolutionary algorithm for constrained MOPs

%------------------------------- Reference --------------------------------
% K. Li, R. Chen, G. Fu, and X. Yao, Two-archive evolutionary algorithm for
% constrained multi-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2018, 23(2): 303-315.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);

            %% Generate random population
            Population = Problem.Initialization();
            CA = UpdateCA([],Population,W);
            DA = UpdateDA(CA,[],Population,W);
            
            CPS = [];
            Metric=[];
            Metric_HV=[];
            igds=[];
            
            objs=[];
            t=1;
            %% Optimization
            while Algorithm.NotTerminated(CA)
                %% mating pool choosing
                % calculate the ratio of non-dominated solutions of CA and DA in Hm
                
                if Changed(Problem,CA)
                    % Save the population before the change

                    Optimums=Problem.Optimums;
                    Optimum_t=cell2mat(Optimums(t,2));
                    IGD_t=IGD(CA,Optimum_t);
                    if isnan(IGD_t)
                        IGD_t=1;
                    end
                    Metric=[Metric;IGD_t];
                    HV_t=HV(CA,Optimum_t);
                    if isnan(HV_t)
                        HV_t=0;
                    end
                    Metric_HV=[Metric_HV,HV_t];
                    t1 = floor((Problem.FE/Problem.N-Problem.preEvolution)/Problem.taut)/Problem.nt;
                    
                    objs1=CA.objs+t1;
%                     objs1=CA.objs+2*t1;
                    objs=[objs;objs1];
                    CPS = [CPS,CA];
                    % React to the change
                    CA = ReinitializationCA(Problem,CA);
                    DA = ReinitializationCA(Problem,DA);
                    t=t+1;
                    
                end
                
                
                Hm = [CA,DA];                         
                [FrontNo,~] = NDSort(Hm.objs,inf);
                FrontNo_C   = FrontNo(1:ceil(length(Hm)/2));
                Nc = size(find(FrontNo_C==1),2);      
                Pc = Nc/length(Hm);
                FrontNo_D = FrontNo(ceil(length(Hm)/2)+1:length(Hm));
                Nd = size(find(FrontNo_D==1),2);      
                Pd = Nd/length(Hm);

                % calculate the proportion of non-dominated solutions in CA
                [FrontNo,~] = NDSort(CA.objs,inf);
                NC = size(find(FrontNo==1),2);         
                PC = NC/length(CA);                     % PC denotes the proportion of non-dominated solutions in CA,it is different from Pc

                %reproduction
                Q = [];
                for i = 1 : size(W,1)
                    if Pc > Pd
                        P1 = MatingSelection(CA); 
                    else
                        P1 = MatingSelection(DA);
                    end
                    pf = rand();
                    if pf < PC
                        P2 = MatingSelection(CA);
                    else
                        P2 = MatingSelection(DA);
                    end
                    MatingPool = [P1,P2];
                    Offspring  = OperatorGAhalf(Problem,MatingPool);
                    Q = [Q,Offspring];
                end

               %% update CA and DA
                CA = UpdateCA(CA,Q,W);
                DA = UpdateDA(CA,DA,Q,W);
                Optimums=Problem.Optimums;
                Optimum_t=cell2mat(Optimums(t,2));
                igd=IGD(CA,Optimum_t);
                igds=[igds,igd];
                if Problem.FE >= Problem.maxFE
                    % Return all populations
                    CPS = [CPS,Population];
                    
                    Optimums=Problem.Optimums;
                    Optimum_t=cell2mat(Optimums(t,2));
                    IGD_t=IGD(CA,Optimum_t);
                    if isnan(IGD_t)
                        IGD_t=1;
                    end
                    Metric=[Metric;IGD_t];
                    HV_t=HV(CA,Optimum_t);
                    if isnan(HV_t)
                        HV_t=0;
                    end
                    Metric_HV=[Metric_HV,HV_t];
                    %文件名如何自定义
                    folder=fullfile('Algorithms\Multi-objective optimization\DCMOEA-CTD\Data',class(Algorithm));
                    [~,~]=mkdir(folder);
                    file   = fullfile(folder,sprintf('%s_',class(Problem)));
                    runNo  = 1;
                    sum_MIGD=0;
                    MIGDs=[];
                    sum_MHV=0;
                    MHVs=[];
                    while exist([file,num2str(runNo),'.mat'],'file') == 2
                        A=load([file,num2str(runNo)]);
%                         MIGD_runNo=load([file,num2str(runNo)],'MIGD');
                        sum_MIGD=sum_MIGD+A.MIGD;
                        sum_MHV=sum_MHV+A.MHV;
                        MIGDs=[MIGDs,A.MIGD];
                        MHVs=[MHVs,A.MHV];
                        runNo = runNo + 1;
                    end
                    %MIGD的均值
                    MIGD=mean(Metric);
                    MIGDs=[MIGDs,MIGD];
                    MIGD_means=(sum_MIGD+MIGD)/runNo;
                    
                    MHV=mean(Metric_HV);
                    MHVs=[MHVs,MHV];
                    MHV_means=(sum_MHV+MHV)/runNo;
                    Std_HV=std(MHVs,1);

                    Std=std(MIGDs,1);
                    save([file,num2str(runNo),'.mat'],'CPS','Metric','MIGD','MIGD_means','MIGDs','Std','Metric_HV','MHV','MHV_means','MHVs','Std_HV','igds');
                end
            end
        end
    end
end