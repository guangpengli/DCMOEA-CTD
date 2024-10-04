function [Population1,Population2,evalues,Fitness1,Fitness2]=Response6(Pop1,Pop2,Problem,CPS,t)
    % type2:目标固定，约束变化
    %% Pop1响应
    S=CalculateSeverity(Problem,Pop1);
    [idxt,Cst]=kmeans(Pop1.decs,S);
    Pop1t_1=CPS(t-1,:);
    [~,Cst_1]=kmeans(Pop1t_1.decs,S);
    distances=pdist2(Cst,Cst_1);
    Population1=Pop1;
    for i=1:S
        cit=Cst(i,:);
        [~,index]=min(distances(i,:));
        cit_1=Cst_1(index,:);
        Deltai_t=cit-cit_1;
        Xi=Population1(idxt==i).decs;
        Xi=Xi+repmat(Deltai_t,size(Xi,1),1);
        Population1(idxt==i)=Problem.Evaluation(Xi);
    end
    Pop2=Problem.Evaluation(Pop2.decs);
    evalues=200;

    [Population1,Fitness1]=EnvironmentalSelection([Population1,Pop2],Problem.N,1);
    [Population2,Fitness2]=EnvironmentalSelection([Population1,Pop2],Problem.N,2);
end