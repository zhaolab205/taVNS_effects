clear
clc
%%
%time_point
load("D:\Fear\analysis\NIRS\data_6sdelay\behaviour_data_s.mat");
load("D:\Fear\analysis\NIRS\data_6sdelay\behaviour_data_t.mat");
P = 0.1;
length_data = size(behaviour_data_s,1);
p_time = zeros(length_data,1);
t_time = zeros(length_data,1);
for i = 1:length_data
    [h,p,ci,t] = ttest2(behaviour_data_s(i,:),behaviour_data_t(i,:));
    p_time(i) = p;
    t_time(i) = t.tstat;
end
place = find(p_time<P);
%%
%data import
load("D:\Fear\analysis\NIRS\data_6sdelay\channel_task_resample_s7.5.mat");
load("D:\Fear\analysis\NIRS\data_6sdelay\channel_task_resample_t7.5.mat");

load("D:\Fear\analysis\NIRS\data_6sdelay\behaviour_data_s.mat");
load("D:\Fear\analysis\NIRS\data_6sdelay\behaviour_data_t.mat");

load('D:\Fear\analysis\NIRS\data_6sdelay\dcc_s.mat');
load('D:\Fear\analysis\NIRS\data_6sdelay\dcc_t.mat');

nchannel = size(channel_task_resample_s,3);
nsub_s = size(channel_task_resample_s,1);
nsub_t = size(channel_task_resample_t,1);
task_length = size(channel_task_resample_s,2);
%%
%channle activation
time_point1 = 14:32;
time_point2 = 237:285;
time_point3 = 296:316;
time_point4 = 1115:1131;
p_act_task = zeros(4,nchannel);
t_act_task = zeros(4,nchannel);
for i =1:nchannel
    [h1,p1,ci1,t1] = ttest2(mean(channel_task_resample_s(:,time_point1,i),2),mean(channel_task_resample_t(:,time_point1,i),2));
    [h2,p2,ci2,t2] = ttest2(mean(channel_task_resample_s(:,time_point2,i),2),mean(channel_task_resample_t(:,time_point2,i),2));
    [h3,p3,ci3,t3] = ttest2(mean(channel_task_resample_s(:,time_point3,i),2),mean(channel_task_resample_t(:,time_point3,i),2));
    [h4,p4,ci4,t4] = ttest2(mean(channel_task_resample_s(:,time_point4,i),2),mean(channel_task_resample_t(:,time_point4,i),2));

    p_act_task(1,i) = p1; t_act_task(1,i) = t1.tstat;
    p_act_task(2,i) = p2; t_act_task(2,i) = t2.tstat;
    p_act_task(3,i) = p3; t_act_task(3,i) = t3.tstat;
    p_act_task(4,i) = p4; t_act_task(4,i) = t4.tstat;

end

mean_diff = @(g1, g2) mean(g1) - mean(g2);
n_permutations = 10000;
for i =1:nchannel
    [p1, perm_stats1] = permutation_test(mean(channel_task_resample_s(:,time_point1,i),2), ...
        mean(channel_task_resample_t(:,time_point1,i),2),n_permutations,mean_diff,'both');
    [p2, perm_stats2] = permutation_test(mean(channel_task_resample_s(:,time_point2,i),2), ...
        mean(channel_task_resample_t(:,time_point2,i),2),n_permutations,mean_diff,'both');
    [p3, perm_stats3] = permutation_test(mean(channel_task_resample_s(:,time_point3,i),2), ...
        mean(channel_task_resample_t(:,time_point3,i),2),n_permutations,mean_diff,'both');
    [p4, perm_stats4] = permutation_test(mean(channel_task_resample_s(:,time_point4,i),2), ...
        mean(channel_task_resample_t(:,time_point4,i),2),n_permutations,mean_diff,'both');
    p_act_task(1,i) = p1; %t_act_task(1,i) = perm_stats1;
    p_act_task(2,i) = p2; %t_act_task(2,i) = perm_stats1;
    p_act_task(3,i) = p3; %t_act_task(3,i) = perm_stats1;
    p_act_task(4,i) = p4; %t_act_task(4,i) = perm_stats1;
end
%relevance to behavioural outcomes
r_act_behaviour_s = zeros(4,nchannel);
r_act_behaviour_t = zeros(4,nchannel);
p_act_behaviour_s = zeros(4,nchannel);
p_act_behaviour_t = zeros(4,nchannel);
for i =1:nchannel
    [r1_s,p1_s] = corr(reshape(channel_task_resample_s(:,time_point1,i),[],1),reshape(behaviour_data_s(time_point1,:)',[],1));
    [r2_s,p2_s] = corr(reshape(channel_task_resample_s(:,time_point2,i),[],1),reshape(behaviour_data_s(time_point2,:)',[],1));
    [r3_s,p3_s] = corr(reshape(channel_task_resample_s(:,time_point3,i),[],1),reshape(behaviour_data_s(time_point3,:)',[],1));
    [r4_s,p4_s] = corr(reshape(channel_task_resample_s(:,time_point4,i),[],1),reshape(behaviour_data_s(time_point4,:)',[],1));
    [r1_t,p1_t] = corr(reshape(channel_task_resample_t(:,time_point1,i),[],1),reshape(behaviour_data_t(time_point1,:)',[],1));
    [r2_t,p2_t] = corr(reshape(channel_task_resample_t(:,time_point2,i),[],1),reshape(behaviour_data_t(time_point2,:)',[],1));
    [r3_t,p3_t] = corr(reshape(channel_task_resample_t(:,time_point3,i),[],1),reshape(behaviour_data_t(time_point3,:)',[],1));
    [r4_t,p4_t] = corr(reshape(channel_task_resample_t(:,time_point4,i),[],1),reshape(behaviour_data_t(time_point4,:)',[],1));
    r_act_behaviour_s(1,i) = r1_s; r_act_behaviour_t(1,i) = r1_t;
    r_act_behaviour_s(2,i) = r2_s; r_act_behaviour_t(2,i) = r2_t;
    r_act_behaviour_s(3,i) = r3_s; r_act_behaviour_t(3,i) = r3_t;
    r_act_behaviour_s(4,i) = r4_s; r_act_behaviour_t(4,i) = r4_t;
    p_act_behaviour_s(1,i) = p1_s; p_act_behaviour_t(1,i) = p1_t;
    p_act_behaviour_s(2,i) = p2_s; p_act_behaviour_t(2,i) = p2_t;
    p_act_behaviour_s(3,i) = p3_s; p_act_behaviour_t(3,i) = p3_t;
    p_act_behaviour_s(4,i) = p4_s; p_act_behaviour_t(4,i) = p4_t;
end

%%
%DCC
dcc_s = zeros(nsub_s,nchannel,nchannel,task_length);
dcc_t = zeros(nsub_t,nchannel,nchannel,task_length);
for i = 1:nsub_s
    sub_data = squeeze(channel_task_resample_s(i,:,:));
    dcc_s(i,:,:,:) = DCC_jj(sub_data, 'whiten', 'simple', 'doverbose','noflat');
end
for i = 1:nsub_t
    sub_data = squeeze(channel_task_resample_t(i,:,:));
    dcc_t(i,:,:,:) = DCC_jj(sub_data, 'whiten', 'simple', 'doverbose','noflat');
end

time_start = [14 237 296 1115];
time_end = [32 285 316 1131];
p_dcc = zeros(4,15,15);
t_dcc = zeros(4,15,15);
for i = 1:nchannel
    for j = 1:nchannel
        for m = 1:length(time_start)
            dcc_channel_s = mean(squeeze(dcc_s(:,i,j,time_start(m):time_end(m))),2);
            dcc_channel_t = mean(squeeze(dcc_t(:,i,j,time_start(m):time_end(m))),2);
            [h,p,ci,t] = ttest2(dcc_channel_s,dcc_channel_t);
            p_dcc(m,i,j) = p;
            t_dcc(m,i,j) = t.tstat;
        end
    end
end
dcc_channel_s = mean(squeeze(dcc_s(:,9,10,time_start(2):time_end(2))),2);
dcc_channel_t = mean(squeeze(dcc_t(:,9,10,time_start(2):time_end(2))),2);
mean_diff = @(g1, g2) mean(g1) - mean(g2);
n_permutations = 10000;
[p, perm_stats] = permutation_test(dcc_channel_s,dcc_channel_t,n_permutations,mean_diff,'both');

%%
clear
clc
%fnirs fmri compare
load('E:\002Fear\analysis\FMRI\roi_time_series_6.mat');
load('E:\002Fear\analysis\NIRS\data_6sdelay\roi_task_300_s_6.mat');
nroi = size(roi_task_300_s_6,3);
length = size(roi_task_300_s_6,2);
r_act_ff = zeros(1,nroi);
p_act_ff = zeros(1,nroi);
fmri_response = zeros(300,nroi);
fnirs_response = zeros(300,nroi);
for roi = 1:nroi
    act_fnirs = mean(squeeze(roi_task_300_s_6(:,:,roi)),1);
    act_fmri = mean(squeeze(roi_time_series_6(:,:,roi)),1);
    act_fnirs_new = zscore(act_fnirs);
    act_fmri_new = zscore(act_fmri);
    fmri_response(:,roi) = act_fmri;
    fnirs_response(:,roi) = act_fnirs;
    %[r,p] = corr(act_fnirs',act_fmri','type','Pearson');
    [r,p] = corr(act_fnirs_new',act_fmri_new','type','Pearson');
    r_act_ff(1,roi) = r;
    p_act_ff(1,roi) = p;
end
P = [.089
.707
.000
.013
.000
.000];
p_fdr = mafdr(P, 'BHFDR', true); 
p_act_ff_fdr = mafdr(p_act_ff, 'BHFDR', true); 
%%
%scr
clear
clc
load('E:\002Fear\analysis\NIRS\data_6sdelay\amp_s.mat');
load('E:\002Fear\analysis\NIRS\data_6sdelay\amp_t.mat');
load('E:\002Fear\analysis\NIRS\data_6sdelay\amp_TTPs.mat');
load('E:\002Fear\analysis\NIRS\data_6sdelay\amp_TTPt.mat');
p_scr = zeros(1,4);
t_scr = zeros(1,4);
for i = 1:4
    amp_data_s = amp_TTPs(:,i);
    amp_data_t = amp_TTPt(:,i);
    clean_amp_s = amp_data_s(~isnan(amp_data_s));
    clean_amp_t = amp_data_t(~isnan(amp_data_t));
    [h,p,ci,t] = ttest2(clean_amp_s,clean_amp_t);
    p_scr(1,i) = p;
    t_scr(1,i) = t.tstat;
end
for i = 1:4
    amp_data_s = amp_s(:,i);
    amp_data_t = amp_t(:,i);
    [h,p,ci,t] = ttest2(amp_data_s,amp_data_t);
    p_scr(1,i) = p;
    t_scr(1,i) = t.tstat;
end