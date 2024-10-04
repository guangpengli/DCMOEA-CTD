classdef DCP1 < PROBLEM
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
            obj.name='DCP1';
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
            g=1+sum((PopDec(:,2:end)-G).^2,2);
            PopObj(:,1)=g.*PopDec(:,1)+G;
            PopObj(:,2)=g.*(1-PopDec(:,1))+G;
        end
        %% Calculate con
        function PopCon=CalCon(obj,PopDec,PopObj)
            c1 = cos(-0.15.*pi).*PopObj(:,2) - sin(-0.15.*pi).*PopObj(:,1) - (2.*sin(5.*pi.*(sin(-0.15.*pi).*PopObj(:,2) + cos(-0.15.*pi).*PopObj(:,1)))).^6;
            Con = -c1;
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            x1=(0:1/(500-1):1)';
            addtion1 = (0.001:0.001:0.4)';
            G = abs(sin(0.5.*pi.*t));
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            for i=1:length(G)
                pf=[];
                P1=x1+G(i);
                P2=1-x1+G(i);
                addtion2 = repmat(min(P1),length(addtion1),1);
                addtion3 = addtion1 + max(P1);
                P1 = [P1;addtion2;addtion3];
                P2 = [P2;addtion3;addtion2];
                c1 = cos(-0.15*pi)*P2 - sin(-0.15*pi)*P1 - (2*sin(5*pi*(sin(-0.15*pi)*P2 + cos(-0.15*pi)*P1))).^6;
                Con =-c1;
                pf(:,1)=P1;
                pf(:,2)=P2;
                pf(Con>0,:)=[];
                pf(NDSort(pf,1)~=1,:)=[];
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
    end
end













