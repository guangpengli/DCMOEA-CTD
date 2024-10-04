classdef DCF9 < PROBLEM
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
            g = 1 + sum((PopDec(:,2:end)-G).^2,2); 
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
            
            W = 0.5*pi*t;                
            c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
            c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
            c1 = c11 .* c12;
            c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
            c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* PopObj(:,1) - PopObj(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
            c3 = c31 .* c32;
            Con = -[c1 c3];
            PopCon=Con;
        end
        function R=GetOptimum(obj,NS)
            V(:,1) = (0:1/NS:1)';
            V(:,2) = 1 - V(:,1) ;
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = sin(0.5*pi*t);
            obj.Optimums = {};
            for i=1:length(t)
                pf=[];  pf2=[];  
                W = 0.5*pi*t(i);
                pf = V;
                c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                c1a = c11 .* c12;
                c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                c3a = c31 .* c32;
                pf(c1a<0 | c3a<0,:) = [];
                if 0.5*t(i)-floor(0.5*t(i)) >= 0.25 && 0.5*t(i)-floor(0.5*t(i)) <= 0.75
                    pf2(:,1) = [zeros(NS+1,1); 1+V(:,1)];
                    pf2(:,2) = [1+V(:,1); zeros(NS+1,1)];
                    c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                    c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                    c1a = c11 .* c12;
                    c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                    c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                    c3a = c31 .* c32;
                    pf2(c1a<0 | c3a<0,:) = [];
                elseif 0.5*t(i)-floor(0.5*t(i)) < 0.25
                    f12 = -(2.2*(1-cos(W))-1.3)/(cos(W)+sin(W));
                    f11 = (2.2*(1-cos(W))-1.3)/(sin(W)-cos(W));
                    f22 = -(2.2*(1-cos(W))-1.8)/(cos(W)+sin(W));
                    f21 = (2.2*(1-cos(W))-1.8)/(sin(W)-cos(W));
                    pf2 = [V.*repmat([f11 f12],size(V,1),1);V.*repmat([f21 f22],size(V,1),1)];
                    c11 = pf2(:,1) + pf2(:,2) -1;
                    pf2(c11<0,:) = [];
                elseif 0.5*t(i)-floor(0.5*t(i)) > 0.75   
                    f12 = -(2.2*(1-cos(W))-3.1)/(cos(W)+sin(W));
                    f11 = (2.2*(1-cos(W))-3.1)/(sin(W)-cos(W));
                    f22 = -(2.2*(1-cos(W))-2.6)/(cos(W)+sin(W));
                    f21 = (2.2*(1-cos(W))-2.6)/(sin(W)-cos(W));
                    pf2 = [V.*repmat([f11 f12],size(V,1),1);V.*repmat([f21 f22],size(V,1),1)];
                    c11 = pf2(:,1) + pf2(:,2) -1;
                    pf2(c11<0,:) = [];
                end
                pf = [pf; pf2];
                pf(NDSort(pf,1)~=1,:) = [];
%                 pf=pf+2*t(i);
                obj.Optimums(i,:)={t(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
    end
end







