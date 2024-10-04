function S=CalculateSeverity(Problem,Pop1)
    RFt=Problem.CalObj(Pop1.decs);
    Ft=Pop1.objs;
    zmin=min(RFt,[],1);
    zmax=max(RFt,[],1);
    
    u=abs(RFt-Ft)./(zmax);
    mu=mean(u);
    delta=(RFt-Ft)./(zmax);
    DOs=mean(abs(delta-mu),'all');
    
    
    RCVt=Problem.CalCon(Pop1.decs,Pop1.objs);
    RCVt=max(0,RCVt);
    CVt=Pop1.cons;
    CVt=max(0,CVt);
    zmin_c=min(RCVt,[],1);
    zmax_c=max(RCVt,[],1);
    
    for i=1:length(zmax_c)
        if zmax_c(i)==0
            zmax_c(i)=1;
        end
    end
    
    
    u_c=abs(RCVt-CVt)./(zmax_c);
    mu_c=mean(u_c);
    delta=(RCVt-CVt)./(zmax_c);
    DCVs=mean(abs(delta-mu_c),'all');
   
    
    
    
    w1=0.5;w2=0.5;
    severity=w1*DOs+w2*DCVs;
    s1=Problem.M;
    s2=2*Problem.M+1;
    S=s1+ceil(severity*(s2-s1));
    if S <=0
        S=3;
    elseif S>5 
        S=5;
    end
end