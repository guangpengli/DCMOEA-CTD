%储存算法A所有测试问题上的MIGD_means
clc;
clear;
DCF_All_MHV_means=[];
DCF_All_MHV_std=[];
for i=1:10
    str='Algorithms\Multi-objective optimization\DCMOEA-CTD\Data10_20\DCMOEACTD2\';
    file=['DCF',num2str(i),'_5'];
    str=[str,file];
    A=load(str);
    
    MHV_means=sprintf('%.4e',A.MHV_means);
    DCF_All_MHV_means=[DCF_All_MHV_means;MHV_means];
%     DCF_All_MHV_means=[DCF_All_MHV_means;A.MHV_means];
    Std_HV=sprintf('%.4e',A.Std_HV);
    DCF_All_MHV_std=[DCF_All_MHV_std;Std_HV];
%     DCF_All_MHV_std=[DCF_All_MHV_std;A.Std_HV];
end


a='(';
a=repmat(a,10,1);
b=')';
b=repmat(b,10,1);
DCF_All_MHV_means=string(DCF_All_MHV_means);
DCF_All_MHV_std=string(DCF_All_MHV_std);
m=DCF_All_MHV_means+a+DCF_All_MHV_std+b;

DCP_All_MHV_means=[];
DCP_All_MHV_std=[];
for i=1:9
    str='Algorithms\Multi-objective optimization\DCMOEA-CTD\Data10_20\DCMOEACTD2\';
    file=['DCP',num2str(i),'_5'];
    str=[str,file];
    A=load(str);
    MHV_means=sprintf('%.4e',A.MHV_means);
    DCP_All_MHV_means=[DCP_All_MHV_means;MHV_means];
%     DCP_All_MHV_means=[DCP_All_MHV_means;A.MHV_means];
    Std_HV=sprintf('%.4e',A.Std_HV);
    DCP_All_MHV_std=[DCP_All_MHV_std;Std_HV];
%     DCP_All_MHV_std=[DCP_All_MHV_std;A.Std_HV];
end
a='(';
a=repmat(a,9,1);
b=')';
b=repmat(b,9,1);
DCP_All_MHV_means=string(DCP_All_MHV_means);
DCP_All_MHV_std=string(DCP_All_MHV_std);
n=DCP_All_MHV_means+a+DCP_All_MHV_std+b;


a=0;
% DCF_All_MIGD_means=[];
% DCF_All_MIGD_std=[];
% for i=1:10
%     str='Algorithms\Multi-objective optimization\DCMOEA-CTD\Data\DCNSGAII\';
%     file=['DCF',num2str(i),'_5'];
%     str=[str,file];
%     A=load(str);
%     
%     DCF_All_MIGD_means=[DCF_All_MIGD_means;A.MIGD_means];
%     
%     DCF_All_MIGD_std=[DCF_All_MIGD_std;A.Std];
% end
% DCP_All_MIGD_means=[];
% DCP_All_MIGD_std=[];
% for i=1:9
%     str='Algorithms\Multi-objective optimization\DCMOEA-CTD\Data\DCNSGAII\';
%     file=['DCP',num2str(i),'_5'];
%     str=[str,file];
%     A=load(str);
%     
%     DCP_All_MIGD_means=[DCP_All_MIGD_means;A.MIGD_means];
%     
%     DCP_All_MIGD_std=[DCP_All_MIGD_std;A.Std];
% end
% a=0;

