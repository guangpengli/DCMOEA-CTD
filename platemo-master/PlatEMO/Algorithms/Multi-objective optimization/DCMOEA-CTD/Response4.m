function [Population1,Population2,evalues,Fitness1,Fitness2]=Response4(Pop1,Pop2,Problem,CPS,t)
    % type2:目标固定，约束变化；type3:约束变化，UPS固定
    %% Pop2响应
    zmin2=min(Pop2.decs,[],1);
    zmax2=max(Pop2.decs,[],1);
    Gauss2=normrnd(0,1,[Problem.N,1]);
    X2=Pop2.decs+Gauss2.*(zmax2-zmin2);
    Population2=Problem.Evaluation(X2);

    
    %% Pop1响应
    S=CalculateSeverity(Problem,Pop1);
    [idxt,Cst]=kmeans(Pop1.decs,S);
    Pop1t_1=CPS(t-1,:);
    [~,Cst_1]=kmeans(Pop1t_1.decs,S);
    distances=pdist2(Cst,Cst_1);
    Population1=Pop1;
    for i=1:S
        cit=Cst(i,:);
        [d,index]=min(distances(i,:));
        cit_1=Cst_1(index,:);
        Deltai_t=cit-cit_1;
        
        Xi=Population1(idxt==i).decs;
%         Gauss=normrnd(0,d,[size(Xi,1),1]);
%         Xi=Xi+repmat(Deltai_t,size(Xi,1),1)+Gauss;
        Xi=Xi+repmat(Deltai_t,size(Xi,1),1);
        Population1(idxt==i)=Problem.Evaluation(Xi);
    end

    %% 选择

    evalues=200;
   [Population1,Fitness1]=EnvironmentalSelection([Population1,Population2],Problem.N,1);
   [Population2,Fitness2]=EnvironmentalSelection([Population1,Population2],Problem.N,2);
end