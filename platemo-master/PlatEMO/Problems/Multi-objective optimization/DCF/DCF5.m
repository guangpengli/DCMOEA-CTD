classdef DCF5 < PROBLEM
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
%             obj.lower    = zeros(1,obj.D);
%             obj.upper    = ones(1,obj.D);
%             obj.encoding = ones(1,obj.D);
            obj.lower    = [0,-ones(1,obj.D-1)];
            obj.upper    = [1, ones(1,obj.D-1)];
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
            G = sin(0.5.*pi.*t);
            
            g = 1+sum((PopDec(:,2:end)- G).^2 +sin((PopDec(:,2:end)- G)*0.5*pi).^2,2);
            PopObj(:,1) = g.*(PopDec(:,1) + 0.2*G*sin(pi*PopDec(:,1)));
            PopObj(:,2) = g.*(1-PopDec(:,1) + 0.2*G*sin(pi*PopDec(:,1)));
            PopObj(PopObj < 1e-18) = 0;
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            c11 = PopObj(:,1) + 2*PopObj(:,2) - 1;
            c12 = PopObj(:,1) + 0.5*PopObj(:,2) - 0.5;
            c1 = c11 .* c12;

            c2 = PopObj(:,1).^2 + PopObj(:,2).^2 - 1.4.^2;
            Con = [-c1 c2];
            PopCon=Con;
        end
        function R=GetOptimum(obj,NS)
            x1 = (0:1/NS:1)';
            pf2=[]; pf2(:,1)=x1;
            pf2(:,2) = [1-2*x1(x1<1/3); 0.5-0.5*x1(x1>=1/3)];
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = sin(0.5*pi*t);
            obj.Optimums = {};
            for i=1:length(t)
                pf=[];  
                pf(:,1) = x1 + 0.2*G(i)*sin(pi*x1);
                pf(:,2) = 1-x1 + 0.2*G(i)*sin(pi*x1);
                pf(pf==min(pf)) = 0;  pf(pf==max(pf)) = 1;
                c11 = pf(:,1) + 2*pf(:,2) - 1;
                c12 = pf(:,1) + 0.5*pf(:,2) - 0.5;
                c1 = c11 .* c12;
                pf(c1<0,:) = [];
                LengthPF = size(pf,1);
                pf = [pf; pf2];
                Select = NDSort(-pf,1) == inf;
                Select(1:LengthPF) = false;
                pf(Select,:) = [];
                pf(NDSort(pf,1)~=1,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={t(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
    end
end







