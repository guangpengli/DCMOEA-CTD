classdef DCP4 < PROBLEM
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
            obj.name='DCP4';
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
            G = 2 * floor(10*abs(mod(t+1,2)-1)+1e-4);
            y = (PopDec(:,2:end)-0.5);
            g = 1 + sum((y.^2 - cos(pi*y) + 1).^2,2);
            
            PopObj(:,1) = g.*(PopDec(:,1) + 0.25*sin(pi*PopDec(:,1)));
            PopObj(:,2) = g.*(1 - PopDec(:,1) + 0.25*sin(pi*PopDec(:,1)));
                
        end
         %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            if obj.FE/obj.N<obj.preEvolution
                t=0;
            else
                t = floor((obj.FE/obj.N-obj.preEvolution)/obj.taut+1)/obj.nt;
            end
            G = 2 * floor(10*abs(mod(t+1,2)-1)+1e-4);
            c11 = (PopObj(:,1)).^2 + (PopObj(:,2)).^2 - (1.5 + 0.4*sin(4*atan(PopObj(:,2)./PopObj(:,1))).^16).^2;
            c12 = (PopObj(:,1)).^2 + (PopObj(:,2)).^2 - (1.3 - 0.45*sin(G*atan(PopObj(:,2)./PopObj(:,1))).^2).^2;
            Con = -c11.*c12;
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = 2 * floor(10*abs(mod(t+1,2)-1)+1e-4);
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            x1 = (0:1/(500-1):1)';
            pf1(:,1) = x1 +  0.25*sin(pi*x1);
            pf1(:,2) = 1- x1 +  0.25*sin(pi*x1);
            
            for i=1:length(G)
                pf=pf1; 
                c11 = (pf(:,1)).^2 + (pf(:,2)).^2 - (1.5 + 0.4*sin(4*atan(pf(:,2)./pf(:,1))).^16).^2;
                c12 = (pf(:,1)).^2 + (pf(:,2)).^2 - (1.3 - 0.45*sin(G(i)*atan(pf(:,2)./pf(:,1))).^2).^2;
                c1 = c11.*c12<0;
                pf(c1,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
       
    end
end
















