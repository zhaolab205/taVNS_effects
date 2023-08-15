function DFC_calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcluate functional connectivity connectome for each subject
% written by
% Qi Liu; Siyu Zhu
% Weihua Zhao: zarazhao@uestc.edu.cn
% last edited Aug 2023

%input--- channel_signal：original signal
%     --- onset_inf     ：onsets information
%     --- channel_seq   ：channel sequence
%                         left IFG: 23 24 25 27;
%                         right IFG: 8 9 10 11;
%                         left DLPFC: 19 20 21 22;
%                         right DLPFC: 12 13 14 15;
%                         medial PFC: 5 6 16 17 18 28 29；
%                         left OFC:  3 26 30 2;
%                         right OFC: (2) 1  4 7；  此处的2为重叠channel

%output--- DFC_cha_inf:82 subjects  30roi*30roi*209win  r_matrix  p_matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load('mean_signal.mat')
load('channel_signal.mat')
load('onset_inf.mat')
load('channel_seq.mat')
[nsub,ncha]=size(channel_signal);
% [nsub,nroi]=size(mean_signal);
nsub=nsub-1;
ncha=30;
nwin = 170;       % window number
step_length = 12; % step length
win_length = 252; %window length

DFC_cha_inf(1,1)={'SubID'};
DFC_cha_inf{1,2}={'r_matrix'};
DFC_cha_inf{1,3}={'p_matrix'};

%% calcluate functional connectivity connectome
for sub_num = 1:nsub
    initial_point = onset_inf{sub_num+1,2}(1,1)-8;  %起始点
    stop_point = onset_inf{sub_num+1,2}(length(onset_inf(sub_num+1,2)),1) ;  %终止点
    DFC_cha_inf{sub_num+1,1}=onset_inf{sub_num+1,1};
    for win_num = 1:nwin
        start_point = initial_point + (win_num-1)*step_length;    %更新初始点
        stop_point = start_point + win_length;
        xroi_num = 0;
       
        for xroi = channel_seq
            xroi_num =  xroi_num +1; 
            yroi_num = 0;
            for yroi = channel_seq
                yroi_num =  yroi_num +1;
                [dfc_r(xroi_num,yroi_num,win_num),dfc_p(xroi_num,yroi_num,win_num)]=corr(channel_signal{sub_num+1,2}(start_point:stop_point,xroi),channel_signal{sub_num+1,2}(start_point:stop_point,yroi));
            end
        end
    end
    DFC_cha_inf{sub_num+1,2} = dfc_r;
    DFC_cha_inf{sub_num+1,3} = dfc_p;
    disp(onset_inf{sub_num+1,1})
end