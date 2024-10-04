classdef DCF6 < PROBLEM
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
            g = 1 + sum((PopDec(:,2:end)-G).^2 ,2); 
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
            G = sin(0.5*pi*t);
            H = 1;
            if G < 0
                H = -1;
            end
            ll  = cos(5*atan((PopObj(:,2)./PopObj(:,1)).^H).^4).^6;
            c12 = 1.1 - (PopObj(:,1)./(0.9 + (0.1+0.7*abs(G))*ll)).^2 - (PopObj(:,2)./(0.9 + (0.8-0.7*abs(G))*ll)).^2 ;
            Con = - c12;
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
                pf=[];  pf = V;
                H = 1;
                if G(i) < 0
                    H = -1;
                end
                ll  = cos(5*atan((pf(:,2)./pf(:,1)).^H).^4).^6;
                Con = 1.1 - (pf(:,1)./(0.9 + (0.1+0.7*abs(G(i)))*ll)).^2 - (pf(:,2)./(0.9 + (0.8-0.7*abs(G(i)))*ll)).^2;
                pf(Con<0,:) = [];
                
                obj.Optimums(i,:)={t(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
    end
end







