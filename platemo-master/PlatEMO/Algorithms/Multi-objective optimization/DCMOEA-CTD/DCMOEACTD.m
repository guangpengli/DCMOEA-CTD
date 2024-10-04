classdef DCMOEACTD < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <dynamic>
% Dynamic DCMOEACTD

    methods
        function main(Algorithm,Problem)
            Algorithm.save = sign(Algorithm.save)*inf;
            Metric_IGD=[];
%             Metric_HV=[];
            CPS=[];
            UPS=[];
            Arch=[];
            % Generate random population
            Population1 = Problem.Initialization();
            Population2 = Population1;
            Fitness1    = CalFitness(Population1,1);
            Fitness2    = CalFitness(Population2,2);
            
            t=1;
            flag=false;
            while Algorithm.NotTerminated(Population1)
                % 动态检测响应
                [Change,type]=Changed(Problem,Population1);
                if Change
                    Arch=EnvironmentalSelection([Arch,Population1,Population2],Problem.N,3);
                    Optimums=Problem.Optimums;
                    Optimum_t=cell2mat(Optimums(t,2));
                    IGD_t=IGD(Arch,Optimum_t);
                    if isnan(IGD_t)
                        IGD_t=1;
                    end
                    Metric_IGD=[Metric_IGD;IGD_t];
%                     HV_t=HV(Arch,Optimum_t);
%                     if isnan(HV_t)
%                         HV_t=0;
%                     end
%                     Metric_HV=[Metric_HV,HV_t];

                    CPS=[CPS;Arch];
                    UPS=[UPS;Population2];
                    if t>=2
                        evalues=0;
                        % 类别检测响应环境
                        if type==1
                            %约束函数固定，UPS不变，多样性策略
                            if flag
                                %种群2收敛到UPF，采用了SPEA2+CDP搜索，将该问题的UPF应用到下一代
                                Population2=Population2t;
                            end
                            [Population1,Population2,evalues,Fitness1,Fitness2]=Response2(Population1,Population2,Problem);
                        elseif type==4
                            %约束函数固定，UPS变化，Pop2差分预测，Pop1多样性策略
                            [Population1,Population2,evalues,Fitness1,Fitness2]=Response3(Population1,Population2,Problem,UPS,t);
                        elseif type==2
                            % type2:目标固定，约束变化,UPS不变，不响应，直接评价
                            Population2=UPS(1,:);
                            [Population1,Population2,evalues,Fitness1,Fitness2]=Response6(Population1,Population2,Problem,CPS,t);     
                        elseif type==3
                            % type3:约束变化，UPS固定
                            if flag
                                Population2=Population2t;
                            end
                            [Population1,Population2,evalues,Fitness1,Fitness2]=Response4(Population1,Population2,Problem,CPS,t);
                        elseif type==5
                            %目标，约束都发生变化
                            [Population1,Population2,evalues,Fitness1,Fitness2]=Response5(Population1,Population2,Problem,UPS,CPS,t);
                        end
                        %为了公平比较，使得每次环境变化都能进化taut代，删除响应阶段的评价次数evalues。
                        Problem.FE=Problem.FE-evalues;
                        
                    else
                        % 多样性响应
                        [Population1,Population2,evalues,Fitness1,Fitness2]=Response1(Population1,Population2,Problem);
                        Problem.FE=Problem.FE-evalues;
                    end
                      
                    Arch=[];
                    Arch=EnvironmentalSelection([Arch,Population1,Population2],Problem.N,3);
                    flag=false;
                    t=t+1;
                end
                % 静态优化
                OldPopulation2=Population2;
                
                MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));

                MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                index=randperm(Problem.N,Problem.N/2);
                Offspring2  = OperatorDE(Problem,Population2(index),Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                
%                 Offspring2=DE_rand_to_best_1(Problem,Population2,Problem.N/2,Fitness2,Population2);
                Arch=EnvironmentalSelection([Arch,Offspring1,Offspring2],Problem.N,3);
                [Population1,Fitness1] = EnvironmentalSelection([Population1,Arch,Offspring1,Offspring2],Problem.N,1);
                if flag==false
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Arch,Offspring1,Offspring2],Problem.N,2);
                else
                    [Population2,Fitness2] = EnvironmentalSelection([Population2,Arch,Offspring1,Offspring2],Problem.N,3);
                end
                  
                if type==2
                    flag=true;
                elseif type==3||type==1
                    if flag==false
                        %计算目标值差
%                         OldObjs=sum(sum(OldPopulation2.objs));
%                         Objs=sum(sum(Population2.objs));
%                         alpha=abs(Objs-OldObjs)/Objs;
                        %计算HV变化
%                         [FrontNo,~] = NDSort(Population2.objs,length(Population2));
%                         NSg=sum(FrontNo==1);
%                         if NSg==length(Population2)
%                             Optimums=Problem.Optimums;
%                             Optimum_t1=cell2mat(Optimums(t,2));
%                             HV_t1=HV(Population2,Optimum_t1);
%                             HV_t=HV(OldPopulation2,Optimum_t1);
%                             alpha=abs(HV_t1-HV_t)/HV_t;
%                             if alpha<=5e-3
%                                 flag=true;
%                                 Population2t=Population2;
%                             end
%                         end

                        Optimums=Problem.Optimums;
                        Optimum_t1=cell2mat(Optimums(t,2));
                        HV_t1=HV(Population2,Optimum_t1);
                        HV_t=HV(OldPopulation2,Optimum_t1);
                        alpha=abs(HV_t1-HV_t)/HV_t;
                        if alpha<=5e-3
                            flag=true;
                            Population2t=Population2;
                        end
                    end
                end

                
               
                if Problem.FE >= Problem.maxFE
                    % Return all populations
                    Population1=EnvironmentalSelection([Arch,Population1,Population2],Problem.N,3);
                    CPS = [CPS;Population1];
                    
                    Optimums=Problem.Optimums;
                    Optimum_t=cell2mat(Optimums(t,2));
                    IGD_t=IGD(Population1,Optimum_t);
                    
                    Metric_IGD=[Metric_IGD;IGD_t];

                    %文件名如何自定义
                    folder=fullfile('Algorithms\Multi-objective optimization\DCMOEA-CTD\Data',class(Algorithm));
%                     folder=fullfile('Algorithms\Multi-objective optimization\9-30',class(Algorithm));
                    [~,~]=mkdir(folder);
                    file   = fullfile(folder,sprintf('%s_',class(Problem)));
                    runNo  = 1;
                    sum_MIGD=0;
                    MIGDs=[];
                    
                    
                    while exist([file,num2str(runNo),'.mat'],'file') == 2
                        A=load([file,num2str(runNo)]);
        %                         MIGD_runNo=load([file,num2str(runNo)],'MIGD');
                        sum_MIGD=sum_MIGD+A.MIGD;
                        
                        MIGDs=[MIGDs,A.MIGD];
                        
                        runNo = runNo + 1;
                    end
                    %MIGD的均值
                    MIGD=mean(Metric_IGD);
                    MIGDs=[MIGDs,MIGD];
                    MIGD_means=(sum_MIGD+MIGD)/runNo;
                    
                    

                    Std=std(MIGDs,1);
                    save([file,num2str(runNo),'.mat'],'Metric_IGD','MIGD','MIGD_means','MIGDs','Std');
                end
           end
        end
        
    end
end

 

   