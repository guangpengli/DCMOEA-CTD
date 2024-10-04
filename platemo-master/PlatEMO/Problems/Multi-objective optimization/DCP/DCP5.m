classdef DCP5 < PROBLEM
% <multi> <real> <large/none> <dynamic>
% Benchmark dynamic MOP proposed by Farina, Deb, and Amato
% taut --- 10 --- Number of generations for static optimization
% nt   --- 10 --- Number of distinct steps
% preEvolution ---60---fdasf
   properties
        taut;       % Number of generations for static optimization
        nt;         % Number of distinct steps
        Optimums;   % Point sets on all Pareto fronts
        preEvolution;
        name;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.name='DCP5';
            [obj.taut,obj.nt,obj.preEvolution] = obj.ParameterSet(10,10,60);
            obj.M = 2;
            if isempty(obj.D); obj.D = 10; end
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
            g = 1 + sum(abs(PopDec(:,2:end) - 0.5*sin(2*pi*PopDec(:,1))),2); 
            PopObj(:,1) = g.*PopDec(:,1);
            PopObj(:,2) = g.*(1 - PopDec(:,1));
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G=0.5.*abs(sin(0.5.*pi.*t));
            c1 = ((0.2+G)*(PopObj(:,1)).^2 + PopObj(:,2) - 2) .* (0.7*(PopObj(:,1)).^2 + PopObj(:,2) - 2.5);
            c2 = (PopObj(:,1)).^2 +(PopObj(:,2)).^2 - (0.6+G ).^2;
            Con = -[c1 c2];
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G=0.5.*abs(sin(0.5.*pi.*t));
%             G = unique(round(G*1e6)/1e6);
            pf1(:,1) = (0:1/(500-1):1)';
            pf1(:,2) = 1- pf1(:,1) ;
            X = UniformPoint(500,2);
            X = X./repmat(sqrt(sum(X.^2,2)),1,2);
            obj.Optimums = {};
            for i=1:length(G)
                pf = pf1;
                c1 = (pf(:,1)).^2 +(pf(:,2)).^2 - (0.6+G(i) ).^2<0;
                pf(c1,:) = [];
                pf = [pf;(0.6+G(i))*X];
                pf(NDSort(-pf,1)~=1,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
    end
end












