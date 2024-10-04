classdef DCF10 < PROBLEM
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
            g = 1 + sum((PopDec(:,2:end)-G).^2 - cos(pi*(PopDec(:,2:end)-G)) +1,2); 
            PopObj(:,1) = g.*PopDec(:,1);
            PopObj(:,2) = g.*sqrt(1-PopDec(:,1).^2);
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G = abs(sin(0.5*pi*t));
            c11 = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.4 + 0.5*abs(G) + (1 - abs(G))*sin(8*atan(PopObj(:,2)./PopObj(:,1))).^12 ).^2;
            c12 = PopObj(:,1).^2 + PopObj(:,2).^2 - (1.4 - (1 - abs(G))*sin(8*atan(PopObj(:,2)./PopObj(:,1))).^12 ).^2;
            Con = -c11.*c12;
            PopCon=Con;
        end
        function R=GetOptimum(obj,NS)
            V(:,1) = (0:1/NS:1)';
            V(:,2) = 1 - V(:,1) ;
            V = V./repmat(sqrt(sum(V.^2,2)),1,2);
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = sin(0.5*pi*t);
            obj.Optimums = {};
            for i=1:length(t)
                pf=[];
                pf = V;
                c11 = pf(:,1).^2 + pf(:,2).^2 - (1.4 + 0.5*abs(G(i)) + (1 - abs(G(i)))*sin(8*atan(pf(:,2)./pf(:,1))).^12 ).^2;
                c12 = pf(:,1).^2 + pf(:,2).^2 - (1.4 - (1 - abs(G(i)))*sin(8*atan(pf(:,2)./pf(:,1))).^12 ).^2;
                Con = c11.*c12;
                pf(Con<0,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={t(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
    end
end







