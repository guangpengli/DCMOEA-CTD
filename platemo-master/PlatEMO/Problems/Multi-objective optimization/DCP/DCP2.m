classdef DCP2 < PROBLEM
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
            obj.name='DCP2';
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
            g=1 + sum((PopDec(:,2:end)- G).^2 + sin((PopDec(:,2:end)- G)*0.5*pi).^2,2);
            PopObj(:,1) = g.*(PopDec(:,1) + 0.25*G*sin(pi*PopDec(:,1)));
            PopObj(:,2) = g.*(1-PopDec(:,1) + 0.25*G*sin(pi*PopDec(:,1)));
            PopObj(PopObj < 1e-18) = 0;
        end
         %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            c1 = (4*PopObj(:,1) + PopObj(:,2) -1).*(0.3*PopObj(:,1) + PopObj(:,2) -0.3);
            c2 = (PopObj(:,1) + PopObj(:,2) -1.3).*(1.85 - PopObj(:,2) - PopObj(:,1) - (0.3*sin(3*pi*(PopObj(:,2)-PopObj(:,1)))).^2 );
            Con = [-c1 c2];
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = sin(0.5.*pi.*t);
            
            x2 = 0 : 0.001 : 2;
            pf2(:,1) = repmat(x2',2001,1);
            for i = 1 : 2001
                pf2(2001*(i-1)+1:2001*i,2) =  0.001 * (i-1);
            end
            cc = (1.85 - pf2(:,2) - pf2(:,1) - (0.3*sin(3*pi*(pf2(:,2)-pf2(:,1)))).^2 );
            pf2(abs(cc) > 2e-3,:) = [];
            x1 = (0:1/(500-1):1)';
            y21 = 1 - 4*x1(x1<0.1891); y22 = 0.3 - 0.3*x1(x1>=0.1891);
            y = [x1 [y21;y22]];
            obj.Optimums={};
            for i=1:length(G)
                pf=[];
                P1 = x1 + 0.25*G(i)*sin(pi*x1);
                P2 = 1-x1 + 0.25*G(i)*sin(pi*x1);
                c1 = (4*P1 + P2 -1).*(0.3*P1 + P2 -0.3)<0;
                c2 = (P1 + P2 -1.3).*(1.85 - P2 - P1 - (0.3*sin(3*pi*(P2-P1))).^2 ) >0;
                pf(:,1) = P1;   pf(:,2) = P2;
                pf(c1 | c2,:) = [];
                pf = [pf;y];
                Select = NDSort(-pf,1)~=1;
                pf(Select,:) = []; [~,MinPF]=min(pf(:,1)); pf(MinPF,1) = 0;
                [~,MinPF]=min(pf(:,2)); pf(MinPF,2) = 0;
                pf = [pf;pf2];
                pf(NDSort(pf,1)~=1,:) = [];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
   end
end













