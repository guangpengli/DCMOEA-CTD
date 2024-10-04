classdef DCNSGAIII < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
% Nondominated sorting genetic algorithm III

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
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
            %% Generate the reference points and random population
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
            CPS = [];
            Metric=[];
            Metric_HV=[];
            t=1;
            igds=[];
            initPops=[];
            objs=[];
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    % Save the population before the change

                    Optimums=Problem.Optimums;
                    Optimum_t=cell2mat(Optimums(t,2));
                    IGD_t=IGD(Population,Optimum_t);
                    if isnan(IGD_t)
                        IGD_t=1;
                    end
                    Metric=[Metric;IGD_t];
                    HV_t=HV(Population,Optimum_t);
                    if isnan(HV_t)
                        HV_t=0;
                    end
                    Metric_HV=[Metric_HV,HV_t];
                    t1 = floor((Problem.FE/Problem.N-Problem.preEvolution)/Problem.taut)/Problem.nt;
                    
                    objs1=Population.objs+t1;
%                     objs1=Population.objs+2*t1;
                    objs=[objs;objs1];
                    CPS = [CPS,Population];
                    % React to the change
                    Population = Reinitialization(Problem,Population,Z);
                    t=t+1;
                    initPops=[initPops;Population];
                end
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                Optimums=Problem.Optimums;
                Optimum_t=cell2mat(Optimums(t,2));
                igd=IGD(Population,Optimum_t);
                igds=[igds,igd];
                
                if Problem.FE >= Problem.maxFE
                    % Return all populations
                    CPS = [CPS,Population];
                    
                    Optimums=Problem.Optimums;
                    Optimum_t=cell2mat(Optimums(t,2));
                    IGD_t=IGD(Population,Optimum_t);
                    if isnan(IGD_t)
                        IGD_t=1;
                    end
                    Metric=[Metric;IGD_t];
                    HV_t=HV(Population,Optimum_t);
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
                    save([file,num2str(runNo),'.mat'],'CPS','Metric','MIGD','MIGD_means','MIGDs','Std','Metric_HV','MHV','MHV_means','MHVs','Std_HV','igds','initPops');
                end
            end
        end
    end
end