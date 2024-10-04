classdef DCP8 < PROBLEM
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
            obj.name='DCP8';
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
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G=sin(0.5.*pi.*t);
            y = (PopDec(:,2:end)-G);
            g = 1 + sum((abs(G)*y.^2 - cos(pi*y) +1).^2,2);
            PopObj(:,1) = g.*PopDec(:,1);
            PopObj(:,2) = g.*(1-PopDec(:,1));
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G=sin(0.5.*pi.*t);
            c11 = (PopObj(:,1)).^1.5+ (PopObj(:,2)).^1.5 - 1.2^1.5;
            c12 = PopObj(:,1).^(0.5) + PopObj(:,2).^(0.5) - (0.95+0.5*abs(G) );
            c21 =  0.8*PopObj(:,1) + PopObj(:,2) - (2.5 + 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1))));
            c22 = (0.93 + abs(G)/3)*PopObj(:,1) + PopObj(:,2) - (2.7+abs(G)/2 + 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1))));
            c1  = c11.*c12;
            c2  = c21.*c22;
            Con = -[c1 c2];
            PopCon=Con;
        end
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G=sin(0.5.*pi.*t);
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            x1 = (0:1/(501-1):1)';
            X = UniformPoint(500,2);
            pf1 = 1.2* X./repmat((sum(X.^1.5,2)).^(2/3),1,2);
            for i=1:length(G)
                pf=[];
                pf(:,1) = x1 ;
                pf(:,2) = 1-x1;
                c1 = pf(:,1).^(0.5) + pf(:,2).^(0.5) - (0.95+0.5*abs(G(i)) );
                pf(c1>0,:) = [];
                pf = [pf;pf1];
                pf(NDSort(pf,1)~=1,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
    end
end
        
        
        
        
        
        
        