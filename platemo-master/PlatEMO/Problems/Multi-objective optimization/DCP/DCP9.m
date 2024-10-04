classdef DCP9 < PROBLEM
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
            obj.name='DCP9';
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
            G=abs(sin(0.5.*pi.*t));
            r = 1 + floor((obj.D-1)*G);
            UnDec = PopDec;
            UnDec(:,r) = [];
            g = 1 + 10*sum((UnDec-G).^2,2); 
            PopObj(:,1) = g.*PopDec(:,r);
            PopObj(:,2) = g.*(1-PopDec(:,r)) ;
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G=abs(sin(0.5.*pi.*t));
            c11 = PopObj(:,1).^2+ PopObj(:,2).^2 - (0.2+G).^2;
            c12 = PopObj(:,1).^(0.75+1.25*G) + PopObj(:,2).^(0.75+1.25*G) - 4.^(0.75+1.25*G);
            c1 = c11.*c12;
            c21 = PopObj(:,1).^2+ PopObj(:,2).^2 - 1.6.^2;
            c22 = 2.1 - (PopObj(:,1)./(1+0.15*cos(6*atan(PopObj(:,2)./PopObj(:,1)).^3).^10)).^2 - (PopObj(:,2)./(1+0.75*cos(6*atan(PopObj(:,2)./PopObj(:,1)).^3).^10)).^2;
            c2 = c21.*c22;
            Con = [c1 c2];
            PopCon=Con;
            
        end
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G=abs(sin(0.5.*pi.*t));
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            x1 = (0:1/(501-1):1)';
            X = UniformPoint(500,2);
            pf1 = X./repmat(sqrt(sum(X.^2,2)),1,2);
            for i=1:length(G)
                pf=[];
                pf(:,1) = x1 ;
                pf(:,2) = 1-x1;
                c1 = (pf(:,1)).^2+ (pf(:,2)).^2 - ((0.2+abs(G(i)))).^2;
                pf(c1<0,:) = [];
                pf = [pf;(0.2+abs(G(i)))*pf1];
                pf(NDSort(-pf,1)~=1,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
       
    end
end











