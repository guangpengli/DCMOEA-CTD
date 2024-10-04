%储存算法A所有测试问题上的MIGD_means
clc;
clear;
DCF_All_MIGD_means=[];
DCF_All_MIGD_std=[];
for i=1:10
    str='Algorithms\Multi-objective optimization\DCMOEA-CTD\variant10_20\DCMOEACTDv3\';
    file=['DCF',num2str(i),'_5'];
    str=[str,file];
    A=load(str);
    
    MIGD_means=sprintf('%.4e',A.MIGD_means);
    DCF_All_MIGD_means=[DCF_All_MIGD_means;MIGD_means];
%     DCF_All_MIGD_means=[DCF_All_MIGD_means;A.MIGD_means];
    Std=sprintf('%.4e',A.Std);
    DCF_All_MIGD_std=[DCF_All_MIGD_std;Std];
%     DCF_All_MIGD_std=[DCF_All_MIGD_std;A.Std_IGD];
end


a='(';
a=repmat(a,10,1);
b=')';
b=repmat(b,10,1);
DCF_All_MIGD_means=string(DCF_All_MIGD_means);
DCF_All_MIGD_std=string(DCF_All_MIGD_std);
m=DCF_All_MIGD_means+a+DCF_All_MIGD_std+b;

DCP_All_MIGD_means=[];
DCP_All_MIGD_std=[];
for i=1:9
    str='Algorithms\Multi-objective optimization\DCMOEA-CTD\variant10_20\DCMOEACTDv3\';
    file=['DCP',num2str(i),'_5'];
    str=[str,file];
    A=load(str);
    MIGD_means=sprintf('%.4e',A.MIGD_means);
    DCP_All_MIGD_means=[DCP_All_MIGD_means;MIGD_means];
%     DCP_All_MIGD_means=[DCP_All_MIGD_means;A.MIGD_means];
    Std=sprintf('%.4e',A.Std);
    DCP_All_MIGD_std=[DCP_All_MIGD_std;Std];
%     DCP_All_MIGD_std=[DCP_All_MIGD_std;A.Std_IGD];
end
a='(';
a=repmat(a,9,1);
b=')';
b=repmat(b,9,1);
DCP_All_MIGD_means=string(DCP_All_MIGD_means);
DCP_All_MIGD_std=string(DCP_All_MIGD_std);
n=DCP_All_MIGD_means+a+DCP_All_MIGD_std+b;


a=0;



















