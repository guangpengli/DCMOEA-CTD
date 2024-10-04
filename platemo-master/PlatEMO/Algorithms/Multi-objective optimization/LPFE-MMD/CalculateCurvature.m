function p=CalculateCurvature(pobjs)
    MMDs=sum(pobjs,2);
    pmin=1;
    MMMDs=MMDs;
    concave=find(MMMDs>pmin);
    MMMDs(concave)=2-MMMDs(concave);
    [~,index]=min(MMMDs);
    p=MMDs(index);
end