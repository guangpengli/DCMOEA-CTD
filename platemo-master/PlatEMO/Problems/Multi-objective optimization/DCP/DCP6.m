classdef DCP6 < PROBLEM
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
            obj.name='DCP6';
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
            g = 1 + 6*sum(PopDec(:,2:end).^2,2); 
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
            G=abs(sin(0.5.*pi.*t)).^0.5;
            c1 = PopObj(:,1) + PopObj(:,2) - (4.5 + 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1)/1.6)));
            c21 = PopObj(:,1) + PopObj(:,2) - (2 - 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1)/1.5)));
            c22 = PopObj(:,1) + PopObj(:,2) - (3.2-G - 0.08*sin(2*pi*(PopObj(:,2)-PopObj(:,1)/1.5)));
            c2 = -c21.*c22;
            c3 = (PopObj(:,1)*cos(-0.25*pi) - PopObj(:,2) *sin(-0.25*pi)).^2 / 1.1.^2 + (PopObj(:,1) *sin(-0.25*pi) +  PopObj(:,2) *cos(-0.25*pi)).^2 /(0.1+G).^2 -(0.1+G).^2; 
            Con = [c1 c2 -c3];
            PopCon=Con;
        end
        %% Generate points on the PF
        function R=GetOptimum(obj,N)
            t = floor(0:(obj.maxFE/obj.N-obj.preEvolution)/obj.taut)/obj.nt;
            G=abs(sin(0.5.*pi.*t)).^0.5;
%             G = unique(round(G*1e6)/1e6);
            obj.Optimums = {};
            x1 = (0:1/(500-1):1)';
            pf1(:,1) = x1 ;
            pf1(:,2) = 1- x1;
            x2 = 0 : 0.001 : 1.5;
            y2(:,1) = repmat(x2',1501,1);
            for i = 1 : 1501
                y2(1501*(i-1)+1:1501*i,2) =  0.001 * (i-1);
            end
            for i=1:length(G)
                pf=pf1;
                pf2 = y2;
                c1 = (pf(:,1)*cos(-0.25*pi)- pf(:,2)*sin(-0.25*pi)).^2/1.1.^2+(pf(:,1)*sin(-0.25*pi)+pf(:,2)*cos(-0.25*pi)).^2/(0.1+G(i)).^2-(0.1+G(i)).^2 <0;
                pf(c1,:) = [];
                c2 = (pf2(:,1)*cos(-0.25*pi)- pf2(:,2)*sin(-0.25*pi)).^2/1.1.^2+(pf2(:,1)*sin(-0.25*pi)+pf2(:,2)*cos(-0.25*pi)).^2/(0.1+G(i)).^2-(0.1+G(i)).^2;
                pf2(c2< -1e-4 | c2> 1e-2,:) = [];
                syms x y
                s=solve((x*cos(-0.25*pi)- y*sin(-0.25*pi)).^2/1.1.^2+(x*sin(-0.25*pi)+y*cos(-0.25*pi)).^2/(0.1+G(i)).^2-(0.1+G(i)).^2,y==1-x,x,y);
                LongX=min(double(s.x));
                LongY=min(double(s.y));
                pf2(pf2(:,1)<LongX | pf2(:,2)<LongY,:)=[];
                pf = [pf;pf2];
                pf(NDSort(pf,1)~=1,:) = [];
                if G(i) == 1
                    pf = [];
                    X = UniformPoint(500,2);
                    pf = (0.1+G(i)).^2 * X./repmat(sqrt(sum(X.^2,2)),1,2);
                end
%                 pf=pf+t(i);
                obj.Optimums(i,:)={G(i),pf};
            end
            R = cat(1,obj.Optimums{:,2});
        end
        
    end
end













