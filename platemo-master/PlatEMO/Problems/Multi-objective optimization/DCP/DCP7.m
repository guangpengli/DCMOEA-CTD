classdef DCP7 < PROBLEM
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
            obj.name='DCP7';
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
            g = 1 + sum((PopDec(:,2:end)-G).^2,2);
            A = 0.02*sin((10-abs(floor(10*G)))*pi*PopDec(:,1));
            PopObj(:,1) = g.*(PopDec(:,1)+ A);
            PopObj(:,2) = g.*(1-PopDec(:,1)+ A);
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G=sin(0.5.*pi.*t);
            c1 = PopObj(:,1) + PopObj(:,2) - sin(5*pi*(PopObj(:,1) - PopObj(:,2) + 1)).^2 - G;
            Con = -c1;
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G=sin(0.5.*pi.*t);
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            x1 = (0:1/(501-1):1)';
            for i=1:length(G)
                pf=[];
                A = 0.02*sin((10-abs(floor(10*G(i))))*pi*x1);
                pf(:,1) = x1 + A;
                pf(:,2) = 1-x1 + A;
                c1 = pf(:,1) + pf(:,2) - sin(5*pi*(pf(:,1) - pf(:,2) + 1)).^2 - G(i);
                pf(c1<0,:) = [];
                pf(NDSort(pf,1)~=1,:)=[];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
    end
end












