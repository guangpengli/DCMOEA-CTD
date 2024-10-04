classdef DCF2 < PROBLEM
% <multi> <real> <large/none> <dynamic>
% Benchmark dynamic MOP proposed by Farina, Deb, and Amato
% taut --- 10 --- Number of generations for static optimization
% nt   --- 10 --- Number of distinct steps
% preEvolution--- --- 第一个环境静态MOEA进化的代数
   properties
        taut;       % Number of generations for static optimization
        nt;         % Number of distinct steps
        Optimums;   % Point sets on all Pareto fronts
        preEvolution;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.taut,obj.nt,obj.preEvolution] = obj.ParameterSet(10,10,60);
            obj.M = 2;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate solutions
        function Population = Evaluation(obj,varargin)
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = obj.CalCon(PopDec,PopObj);
            % Attach the current number of function evaluations to solutions
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);
            obj.FE     = obj.FE + length(Population);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G = abs(sin(0.5*pi*t));
            r = 1 + floor((obj.D-1)*G);
            UnDec = PopDec;
            UnDec(:,r) = [];
            g = 1 + sum((UnDec-G).^2,2); 
            PopObj(:,1) = g.*PopDec(:,r);
            PopObj(:,2) = g.*(1-PopDec(:,r).^2);
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G = abs(sin(0.5*pi*t));
            c1 = cos(-0.15*pi)*PopObj(:,2) - sin(-0.15*pi)*PopObj(:,1) - (2*sin(4*pi*(sin(-0.15*pi)*PopObj(:,2) + cos(-0.15*pi)*PopObj(:,1)))).^6;                
            Con = c1;
            PopCon=Con;
        end
        function R=GetOptimum(obj,NS)
            x1 = (0:1/NS:1)';
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = abs(sin(0.5.*pi.*t));
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            for i=1:length(t)
                pf=[];
                P1 = x1;
                P2 = 1-x1.^2;
                Con = cos(-0.15*pi)*P2 - sin(-0.15*pi)*P1 - (2*sin(4*pi*(sin(-0.15*pi)*P2 + cos(-0.15*pi)*P1))).^6;
                pf(:,1) = P1;
                pf(:,2) = P2;
                pf(Con>0,:) = [];
%                 pf=pf+t(i);
%                 R(i) = struct('PF',pf);
                obj.Optimums(i,:)={t(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        %% Display a population in the objective space
        
%         function DrawObj(obj,Population)
%             if obj.FE/obj.N<obj.preEvolution
%                 t=0;
%             else
%                 t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
%             end
%             tempStream = RandStream('mlfg6331_64','Seed',2);
%             color = rand(tempStream,1,3);
%             Draw(Population.objs,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
%             opt=cell2mat(obj.Optimums(:,1));
%             index=find(opt==t);
%             optimums=cell2mat(obj.Optimums(index,2));
%             Draw(optimums,'*','MarkerSize',1,'Color',color);
%                 
% %             t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
% %             for i=1:length(t)
% %                 tempStream = RandStream('mlfg6331_64','Seed',2);
% %                 color = rand(tempStream,1,3);
% %                 Draw(Population.objs,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
% %                 optimums=cell2mat(obj.Optimums(t(i),2));
% %                 Draw(optimums,'*','MarkerSize',1,'Color',color);
% %             
% %             end
%             
%             
%         end
    end
end
