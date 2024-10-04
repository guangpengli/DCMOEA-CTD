classdef DCP3 < PROBLEM
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
            obj.name='DCP3';
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
            g = 1 + sum(sqrt(abs(PopDec(:,2:end)-G)),2);
            PopObj(:,1) = g.*(PopDec(:,1) )+ G^2;
            PopObj(:,2) = g.*(1-PopDec(:,1) )+ G^2;
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            c1 = ( (PopObj(:,1)-1).^2 + (PopObj(:,2)-0.2).^2 - 0.3^2 ).* ( (PopObj(:,2)-1).^2 + (PopObj(:,1)-0.2).^2 - 0.3^2 );
            c2 = PopObj(:,1).^2 + PopObj(:,2).^2 - 4.^2;
            c3 = (PopObj(:,1).^2 + PopObj(:,2).^2 - (3.1 + 0.2*sin(4*atan(PopObj(:,2)./PopObj(:,1))).^2).^2) .* ((PopObj(:,1)).^2 + (PopObj(:,2)).^2 - (2.3).^2);
            Con = [-c1 c2 -c3];
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G = abs(sin(0.5.*pi.*t));
            
            x1 = (0:1/(500-1):1)';
            x2 = 0 : 0.001 : 2;
            pf2(:,1) = repmat(x2',2001,1);
            for i = 1 : 2001
                pf2(2001*(i-1)+1:2001*i,2) =  0.001 * (i-1);
            end
            cc = ( (pf2(:,1)-1).^2 + (pf2(:,2)-0.2).^2 - (0.3 ).^2 ).* ( (pf2(:,2)-1).^2 + (pf2(:,1)-0.2).^2 - (0.3 ).^2 );
            pf2(abs(cc) > 1e-3 ,:) = [];
            pf2(pf2(:,1)<1 & pf2(:,2)<1,:) = [];
            obj.Optimums = {};
            for i=1:length(G)
                pf=[];
                P1 = x1 + G(i)^2;
                P2 = 1-x1 + G(i)^2;
                pf(:,1) = P1;
                pf(:,2) = P2;
                c1 = ( (pf(:,1)-1).^2 + (pf(:,2)-0.2).^2 - (0.3 ).^2 ).* ( (pf(:,2)-1).^2 + (pf(:,1)-0.2).^2 - (0.3 ).^2 ) < 0;
                pf(c1,:) = [];
                if size(pf,1) < length(x1)
                    pf22 =pf2;
                    pf22(pf22(:,1)<G(i)^2 | pf22(:,2)<G(i)^2,:) = [];
                    De = max(pf(:,1));
                    pf = [pf ;pf22];
                    if 2*De >= 2
                        Select = NDSort(-pf,1)~=1;
                        pf(Select,:) = [];
                    end
                    Select = NDSort(pf,1)~=1;
                    pf(Select,:) = [];
                end
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
    end
end













