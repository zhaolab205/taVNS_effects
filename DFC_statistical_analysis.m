function DFC_statistical_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical Analysis of DFC
% written by
% Qi Liu; Siyu Zhu
% Weihua Zhao: zarazhao@uestc.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('DFC_unfold_inf.mat')
load('IDX_4')
load('DFC_cha_inf')
load('treat1_ID.mat')
load('treat2_ID.mat')
load('channel_seq.mat')
load('gs_v_go.mat')% change behavioral indice
behavioral_index = gs_v_go;% change behavioral indices

nclu=4;
nsub=82;
nwin=170;

DFC_statistic_cal{1,1}={'cluster'};
DFC_statistic_cal{1,2}={'DFC_cluster_mean'};
DFC_statistic_cal{1,3}={'ttest_h'}; % treatment differences
DFC_statistic_cal{1,4}={'ttest_p'};
DFC_statistic_cal{1,5}={'ttest_h_matrix'};
DFC_statistic_cal{1,25}={'ttest_data'};

DFC_statistic_cal{1,6}={'corr_r'}; % correlation results including all subjects
DFC_statistic_cal{1,7}={'corr_p'};
DFC_statistic_cal{1,8}={'corr_h'};
DFC_statistic_cal{1,9}={'corr_h_matrix'};
DFC_statistic_cal{1,10}={'data_corr_all'};

DFC_statistic_cal{1,11}={'corr_r1'}; % treatment1
DFC_statistic_cal{1,12}={'corr_p1'};
DFC_statistic_cal{1,13}={'corr_h1'};
DFC_statistic_cal{1,14}={'corr_h1_matrix'};
DFC_statistic_cal{1,15}={'data_corr1_all'};

DFC_statistic_cal{1,16}={'corr_r2'};% treatment2
DFC_statistic_cal{1,17}={'corr_p2'};
DFC_statistic_cal{1,18}={'corr_h2'};
DFC_statistic_cal{1,19}={'corr_h2_matrix'};
DFC_statistic_cal{1,20}={'data_corr2_all'};
% FC which correlated with behavioral indices and different in treatment
DFC_statistic_cal{1,21}={'overlap_sig_r&t'};
DFC_statistic_cal{1,22}={'overlap_r&t_fc'};
DFC_statistic_cal{1,23}={'overlap_r&t_roi'};
DFC_statistic_cal{1,24}={'overlap_coordination'};
% record FCs correlated with behavioral indices
DFC_fa_significant_cha{1,1}={'cluster'};
DFC_fa_significant_cha{1,2}={'Tol'};
DFC_fa_significant_cha{1,3}={'treatment1'};
DFC_fa_significant_cha{1,4}={'treatment2'};
DFC_statistic_cal{1,26}={'Fisher_z'};
DFC_statistic_cal{1,27}={'Fisher_z_h'};
DFC_statistic_cal{1,28}={'Fisher_z_data'};
DFC_statistic_cal{1,29}={'Fisher_z_roi'};
nper = 5000;
for clu_num = 1:nclu
    DFC_statistic_cal{clu_num+1,1}=clu_num;
    DFC_fa_significant_cha{clu_num+1,1}=clu_num;
    clear sub_cluster_mean
    sub_cluster_mean =zeros(nsub,435);
    for sub_num = 1:nsub
        start_win = ( sub_num - 1) * nwin + 1;
        stop_win = sub_num * nwin;
        cluster_mean = zeros(1,435); 
        cluster_num = 0;
        for win_num = start_win : stop_win
            if IDX4(win_num)==clu_num
                cluster_mean = cluster_mean +DFC_unfold_inf(win_num,2:436);
                cluster_num=cluster_num+1;
            end
        end
        if cluster_num~=0
            cluster_mean = cluster_mean/cluster_num;
        end
        sub_cluster_mean(sub_num,1:435) = cluster_mean;
    end
    DFC_statistic_cal{clu_num+1,1}=clu_num;
    DFC_statistic_cal{clu_num+1,2}=sub_cluster_mean;

    ttest_p=zeros(1,435);
    ttest_h=zeros(1,435);
    corr_r=zeros(1,435);
    corr_p=zeros(1,435);
    corr_h=zeros(1,435);
    corr_r1=zeros(1,435);
    corr_p1=zeros(1,435);
    corr_h1=zeros(1,435);
    corr_r2=zeros(1,435);
    corr_p2=zeros(1,435);
    corr_h2=zeros(1,435);
    fisher_z = zeros(1,435);
    fisher_z_h = zeros(1,435);
    clear data_corr_all data_corr1_all data_corr2_all total_signif total_signif1 total_signif2 total_signif_overlap overlap_coordination
    data_corr_all_num=0;
    corr_ttest_fc_num=0;
    data_corr1_all_num=0;
    data_corr2_all_num=0;
    data_ttest_fc_num=0;
    data_corr1_all = zeros(41,1);
    data_corr2_all = zeros(41,1);
    data_corr_all = zeros(82,1);
    corr_ttest_fc = zeros(82,1);
    total_signif = zeros(2,1);
    total_signif1 = zeros(2,1);
    total_signif2 = zeros(2,1);
    data_ttest_fc=zeros(82,1);
    total_signif_overlap = zeros(2,1);
    overlap_coordination = zeros(2,1);
    for fc_num = 1:435
       disp(fc_num)
       [ttest_h(fc_num),ttest_p(fc_num)]=ttest2(sub_cluster_mean(treat1_ID,fc_num),sub_cluster_mean(treat2_ID,fc_num));
       ttest_p(fc_num) = permuation_test(sub_cluster_mean(:,fc_num),41,41,10000,ttest_p(fc_num));
       [corr_r(fc_num),corr_p(fc_num)]=corr(sub_cluster_mean(:,fc_num),behavioral_index);
       corr_r(fc_num)= permuation_corr(sub_cluster_mean(:,fc_num),behavioral_index,10000,corr_r(fc_num));
       [corr_r1(fc_num),corr_p1(fc_num)]=corr(sub_cluster_mean(treat1_ID,fc_num),behavioral_index(treat1_ID));
       [corr_r2(fc_num),corr_p2(fc_num)]=corr(sub_cluster_mean(treat2_ID,fc_num),behavioral_index(treat2_ID));
       if ttest_h(fc_num)==1
           data_ttest_fc(1:82,data_ttest_fc_num+1) =  sub_cluster_mean(:,fc_num);
           data_ttest_fc_num=data_ttest_fc_num+1;
       end
       if corr_p(fc_num) <= 0.05
           corr_h(fc_num)=1;
           data_corr_all(1:82,data_corr_all_num+1)=sub_cluster_mean(:,fc_num);
           data_corr_all_num = data_corr_all_num+1;
       end
       if corr_p(fc_num) <= 0.05 & ttest_h(fc_num) ==1
           corr_ttest_fc(1:82,corr_ttest_fc_num+1) =  sub_cluster_mean(:,fc_num);
           corr_ttest_fc_num=corr_ttest_fc_num+1;
       end
       if corr_p1(fc_num)<=0.05
           corr_h1(fc_num)=1;
           data_corr1_all(1:41,data_corr1_all_num+1)=sub_cluster_mean(treat1_ID,fc_num);
           data_corr1_all_num = data_corr1_all_num+1;
       end
       if corr_p2(fc_num)<=0.05
           corr_h2(fc_num)=1;
           data_corr2_all(1:41,data_corr2_all_num+1)=sub_cluster_mean(treat2_ID,fc_num);
           data_corr2_all_num = data_corr2_all_num+1;
       end
    end
    
    DFC_statistic_cal{clu_num+1,3}=ttest_h;
    DFC_statistic_cal{clu_num+1,4}=ttest_p;
    
    DFC_statistic_cal{clu_num+1,6}=corr_r;
    DFC_statistic_cal{clu_num+1,7}=corr_p;
    DFC_statistic_cal{clu_num+1,8}=corr_h;
    DFC_statistic_cal{clu_num+1,10}=data_corr_all;
    
    DFC_statistic_cal{clu_num+1,11}=corr_r1;
    DFC_statistic_cal{clu_num+1,12}=corr_p1;
    DFC_statistic_cal{clu_num+1,13}=corr_h1;
    DFC_statistic_cal{clu_num+1,15}=data_corr1_all;
    
    DFC_statistic_cal{clu_num+1,16}=corr_r2;
    DFC_statistic_cal{clu_num+1,17}=corr_p2;
    DFC_statistic_cal{clu_num+1,18}=corr_h2;  
    DFC_statistic_cal{clu_num+1,20}=data_corr2_all;
    
    DFC_statistic_cal{clu_num+1,22}=corr_ttest_fc;
    DFC_statistic_cal{clu_num+1,25}=data_ttest_fc;
    
    DFC_statistic_cal{clu_num+1,26}=(log((corr_r1+1)./(1-corr_r1))/2-log((corr_r2+1)./(1-corr_r2))/2)/(sqrt(1/38+1/38));
    z_high_index_pos = find(DFC_statistic_cal{clu_num+1,26}>1.96);
    z_high_index_neg = find(DFC_statistic_cal{clu_num+1,26}<-1.96);
    fisher_z_h(1,z_high_index_pos)=1;
    fisher_z_h(1,z_high_index_neg)=-1;
    DFC_statistic_cal{clu_num+1,27}=line_to_matrix(fisher_z_h,30,0);
    DFC_statistic_cal{clu_num+1,28}=sub_cluster_mean(:,find(fisher_z_h~=0));
    

    
    
    ttest_h_matrix = zeros(30,30);
    ttest_h_num = 0;
    for xroi = 1:29
        for yroi = xroi+1:30
            ttest_h_matrix(xroi,yroi) = ttest_h(ttest_h_num+1);
            ttest_h_matrix(yroi,xroi) = ttest_h(ttest_h_num+1);
            ttest_h_num=ttest_h_num+1;
        end
    end
    DFC_statistic_cal{clu_num+1,5}=ttest_h_matrix;
    

    corr_h_matrix = zeros(30,30);
    corr_h_num = 0;
    signif_num=0;
    for xroi = 1:29
        for yroi = xroi+1:30
            corr_h_matrix(xroi,yroi) = corr_h(corr_h_num+1);
            corr_h_matrix(yroi,xroi) = corr_h(corr_h_num+1);
            corr_h_num=corr_h_num+1;
            if corr_h(corr_h_num)==1
                signif_num = signif_num+1;
                if xroi==1||xroi==2||xroi==3||xroi==4
                    total_signif(1,signif_num)=1;
                end 
                if xroi==5||xroi==6||xroi==7||xroi==8
                    total_signif(1,signif_num)=2;
                end 
                if xroi==9||xroi==10||xroi==11||xroi==12
                    total_signif(1,signif_num)=3;
                end 
                if xroi==13||xroi==14||xroi==15||xroi==16
                    total_signif(1,signif_num)=4;
                end 
                if xroi==17||xroi==18||xroi==19||xroi==20||xroi==21||xroi==22||xroi==23
                    total_signif(1,signif_num)=5;
                end 
                if xroi==24||xroi==25||xroi==26
                    total_signif(1,signif_num)=6;
                end 
                if xroi==27
                    total_signif(1,signif_num)=67;
                end 
                if xroi==28||xroi==29||xroi==30
                    total_signif(1,signif_num)=7;
                end 

                
                if yroi==1||yroi==2||yroi==3||yroi==4
                    total_signif(2,signif_num)=1;
                end 
                if yroi==5||yroi==6||yroi==7||yroi==8
                    total_signif(2,signif_num)=2;
                end 
                if yroi==9||yroi==10||yroi==11||yroi==12
                    total_signif(2,signif_num)=3;
                end 
                if yroi==13||yroi==14||yroi==15||yroi==16
                    total_signif(2,signif_num)=4;
                end 
                if yroi==17||yroi==18||yroi==19||yroi==20||yroi==21||yroi==22||yroi==23
                    total_signif(2,signif_num)=5;
                end 
                if yroi==24||yroi==25||yroi==26
                    total_signif(2,signif_num)=6;
                end 
                if yroi==27
                    total_signif(2,signif_num)=67;
                end 
                if yroi==28||yroi==29||yroi==30
                    total_signif(2,signif_num)=7;
                end 

            end
        end
    end
    DFC_statistic_cal{clu_num+1,9}=corr_h_matrix;
    DFC_fa_significant_cha{clu_num+1,2} = total_signif;
    

    corr_h1_matrix = zeros(30,30);
    corr_h1_num = 0;
    signif_num1=0;
    for xroi = 1:29
        for yroi = xroi+1:30
            corr_h1_matrix(xroi,yroi) = corr_h1(corr_h1_num+1);
            corr_h1_matrix(yroi,xroi) = corr_h1(corr_h1_num+1);
            corr_h1_num=corr_h1_num+1;
            if corr_h1(corr_h1_num)==1
                signif_num1 = signif_num1+1;
                if xroi==1||xroi==2||xroi==3||xroi==4
                    total_signif1(1,signif_num1)=1;
                end 
                if xroi==5||xroi==6||xroi==7||xroi==8
                    total_signif1(1,signif_num1)=2;
                end 
                if xroi==9||xroi==10||xroi==11||xroi==12
                    total_signif1(1,signif_num1)=3;
                end 
                if xroi==13||xroi==14||xroi==15||xroi==16
                    total_signif1(1,signif_num1)=4;
                end 
                if xroi==17||xroi==18||xroi==19||xroi==20||xroi==21||xroi==22||xroi==23
                    total_signif1(1,signif_num1)=5;
                end 
                if xroi==24||xroi==25||xroi==26
                    total_signif1(1,signif_num1)=6;
                end 
                if xroi==27
                    total_signif1(1,signif_num1)=67;
                end 
                if xroi==28||xroi==29||xroi==30
                    total_signif1(1,signif_num1)=7;
                end 

                
                if yroi==1||yroi==2||yroi==3||yroi==4
                    total_signif1(2,signif_num1)=1;
                end 
                if yroi==5||yroi==6||yroi==7||yroi==8
                    total_signif1(2,signif_num1)=2;
                end 
                if yroi==9||yroi==10||yroi==11||yroi==12
                    total_signif1(2,signif_num1)=3;
                end 
                if yroi==13||yroi==14||yroi==15||yroi==16
                    total_signif1(2,signif_num1)=4;
                end 
                if yroi==17||yroi==18||yroi==19||yroi==20||yroi==21||yroi==22||yroi==23
                    total_signif1(2,signif_num1)=5;
                end 
                if yroi==24||yroi==25||yroi==26
                    total_signif1(2,signif_num1)=6;
                end 
                if yroi==27
                    total_signif1(2,signif_num1)=67;
                end 
                if yroi==28||yroi==29||yroi==30
                    total_signif1(2,signif_num1)=7;
                end 
                
                
            end
        end
    end
    DFC_statistic_cal{clu_num+1,14}=corr_h1_matrix;
    DFC_fa_significant_cha{clu_num+1,3} = total_signif1;
    

    corr_h2_matrix = zeros(30,30);
    corr_h2_num = 0;
    signif_num2=0;
    for xroi = 1:29
        for yroi = xroi+1:30
            corr_h2_matrix(xroi,yroi) = corr_h2(corr_h2_num+1);
            corr_h2_matrix(yroi,xroi) = corr_h2(corr_h2_num+1);
            corr_h2_num=corr_h2_num+1;
            if corr_h2(corr_h2_num)==1
                signif_num2 = signif_num2+1;
                if xroi==1||xroi==2||xroi==3||xroi==4
                    total_signif2(1,signif_num2)=1;
                end
                if xroi==5||xroi==6||xroi==7||xroi==8
                    total_signif2(1,signif_num2)=2;
                end 
                if xroi==9||xroi==10||xroi==11||xroi==12
                    total_signif2(1,signif_num2)=3;
                end 
                if xroi==13||xroi==14||xroi==15||xroi==16
                    total_signif2(1,signif_num2)=4;
                end 
                if xroi==17||xroi==18||xroi==19||xroi==20||xroi==21||xroi==22||xroi==23
                    total_signif2(1,signif_num2)=5;
                end 
                if xroi==24||xroi==25||xroi==26
                    total_signif2(1,signif_num2)=6;
                end 
                if xroi==27
                    total_signif2(1,signif_num2)=67;
                end 
                if xroi==28||xroi==29||xroi==30
                    total_signif2(1,signif_num2)=7;
                end 

                
                if yroi==1||yroi==2||yroi==3||yroi==4
                    total_signif2(2,signif_num2)=1;
                end 
                if yroi==5||yroi==6||yroi==7||yroi==8
                    total_signif2(2,signif_num2)=2;
                end 
                if yroi==9||yroi==10||yroi==11||yroi==12
                    total_signif2(2,signif_num2)=3;
                end 
                if yroi==13||yroi==14||yroi==15||yroi==16
                    total_signif2(2,signif_num2)=4;
                end 
                if yroi==17||yroi==18||yroi==19||yroi==20||yroi==21||yroi==22||yroi==23
                    total_signif2(2,signif_num2)=5;
                end 
                if yroi==24||yroi==25||yroi==26
                    total_signif2(2,signif_num2)=6;
                end 
                if yroi==27
                    total_signif2(2,signif_num2)=67;
                end 
                if yroi==28||yroi==29||yroi==30
                    total_signif2(2,signif_num2)=7;
                end 
                
            end
        end
    end
    
    DFC_statistic_cal{clu_num+1,19}=corr_h2_matrix;
    DFC_fa_significant_cha{clu_num+1,4} = total_signif2;
    
    overlap_sig_r_t = DFC_statistic_cal{clu_num+1,5} + DFC_statistic_cal{clu_num+1,9};
    DFC_statistic_cal{clu_num+1,21} = overlap_sig_r_t;

    signif_overlap=0;
    for xroi = 1:29
        for yroi = xroi+1:30
            if  DFC_statistic_cal{clu_num+1,21}(xroi,yroi)==2
                signif_overlap = signif_overlap+1;
                overlap_coordination(1,signif_overlap) = xroi;
                overlap_coordination(2,signif_overlap) = yroi;
                if xroi==1||xroi==2||xroi==3||xroi==4                   
                    total_signif_overlap(1,signif_overlap)=1;
                end
                if xroi==5||xroi==6||xroi==7||xroi==8
                    total_signif_overlap(1,signif_overlap)=2;
                end 
                if xroi==9||xroi==10||xroi==11||xroi==12
                    total_signif_overlap(1,signif_overlap)=3;
                end 
                if xroi==13||xroi==14||xroi==15||xroi==16
                    total_signif_overlap(1,signif_overlap)=4;
                end 
                if xroi==17||xroi==18||xroi==19||xroi==20||xroi==21||xroi==22||xroi==23
                    total_signif_overlap(1,signif_overlap)=5;
                end 
                if xroi==24||xroi==25||xroi==26
                    total_signif_overlap(1,signif_overlap)=6;
                end 
                if xroi==27
                    total_signif_overlap(1,signif_overlap)=67;
                end 
                if xroi==28||xroi==29||xroi==30
                    total_signif_overlap(1,signif_overlap)=7;
                end 

                
                if yroi==1||yroi==2||yroi==3||yroi==4
                    total_signif_overlap(2,signif_overlap)=1;
                end 
                if yroi==5||yroi==6||yroi==7||yroi==8
                    total_signif_overlap(2,signif_overlap)=2;
                end 
                if yroi==9||yroi==10||yroi==11||yroi==12
                    total_signif_overlap(2,signif_overlap)=3;
                end 
                if yroi==13||yroi==14||yroi==15||yroi==16
                    total_signif_overlap(2,signif_overlap)=4;
                end 
                if yroi==17||yroi==18||yroi==19||yroi==20||yroi==21||yroi==22||yroi==23
                    total_signif_overlap(2,signif_overlap)=5;
                end 
                if yroi==24||yroi==25||yroi==26
                    total_signif_overlap(2,signif_overlap)=6;
                end 
                if yroi==27
                    total_signif_overlap(2,signif_overlap)=67;
                end 
                if yroi==28||yroi==29||yroi==30
                    total_signif_overlap(2,signif_overlap)=7;
                end 
                
            end
        end
    end
    
    DFC_statistic_cal{clu_num+1,23} = total_signif_overlap;
    DFC_statistic_cal{clu_num+1,24} = overlap_coordination;
    

    signif_overlap=0;
    fisherz_roi = zeros(8,1);
    fisher_z_fc_num = 0;
    for xroi = 1:29
        for yroi = xroi+1:30
            fisher_z_fc_num = fisher_z_fc_num +1;
            if  DFC_statistic_cal{clu_num+1,27}(xroi,yroi)==1 || DFC_statistic_cal{clu_num+1,27}(xroi,yroi) == -1
                signif_overlap = signif_overlap+1;
                

                per_z = zeros(nper,1);
                for per_num = 1 : nper
                    per_treat1_ID = treat1_ID(randperm(41));
                    per_treat2_ID = treat2_ID(randperm(41));
                    [per_corr_r1,per_corr_p1]=corr(DFC_statistic_cal{clu_num+1,28}(treat1_ID,signif_overlap),behavioral_index(per_treat1_ID));
                    [per_corr_r2,per_corr_p2]=corr(DFC_statistic_cal{clu_num+1,28}(treat2_ID,signif_overlap),behavioral_index(per_treat2_ID));
                    per_z(per_num,1) = (log((1+per_corr_r1)/(1-per_corr_r1))/2-log((1+per_corr_r2)/(1-per_corr_r2))/2)/sqrt(1/38+1/38);
                    
                end
                fisherz_roi(8,signif_overlap)=min(length(find(per_z>DFC_statistic_cal{clu_num+1,26}(fisher_z_fc_num))),length(find(per_z<DFC_statistic_cal{clu_num+1,26}(fisher_z_fc_num))));
                
                fisherz_roi(1,signif_overlap) = channel_seq(xroi);
                fisherz_roi(2,signif_overlap) = channel_seq(yroi);
                fisherz_roi(5,signif_overlap) = DFC_statistic_cal{clu_num+1,26}(fisher_z_fc_num);
                fisherz_roi(6,signif_overlap) = corr_r1(fisher_z_fc_num);
                fisherz_roi(7,signif_overlap) = corr_r2(fisher_z_fc_num);
                if xroi==1||xroi==2||xroi==3||xroi==4                   
                    fisherz_roi(3,signif_overlap)=1;
                end
                if xroi==5||xroi==6||xroi==7||xroi==8
                    fisherz_roi(3,signif_overlap)=2;
                end 
                if xroi==9||xroi==10||xroi==11||xroi==12
                    fisherz_roi(3,signif_overlap)=3;
                end 
                if xroi==13||xroi==14||xroi==15||xroi==16
                    fisherz_roi(3,signif_overlap)=4;
                end 
                if xroi==17||xroi==18||xroi==19||xroi==20||xroi==21||xroi==22||xroi==23
                    fisherz_roi(3,signif_overlap)=5;
                end 
                if xroi==24||xroi==25||xroi==26
                    fisherz_roi(3,signif_overlap)=6;
                end 
                if xroi==27
                    fisherz_roi(3,signif_overlap)=67;
                end 
                if xroi==28||xroi==29||xroi==30
                    fisherz_roi(3,signif_overlap)=7;
                end 

                
                if yroi==1||yroi==2||yroi==3||yroi==4
                    fisherz_roi(4,signif_overlap)=1;
                end 
                if yroi==5||yroi==6||yroi==7||yroi==8
                    fisherz_roi(4,signif_overlap)=2;
                end 
                if yroi==9||yroi==10||yroi==11||yroi==12
                    fisherz_roi(4,signif_overlap)=3;
                end 
                if yroi==13||yroi==14||yroi==15||yroi==16
                    fisherz_roi(4,signif_overlap)=4;
                end 
                if yroi==17||yroi==18||yroi==19||yroi==20||yroi==21||yroi==22||yroi==23
                    fisherz_roi(4,signif_overlap)=5;
                end 
                if yroi==24||yroi==25||yroi==26
                    fisherz_roi(4,signif_overlap)=6;
                end 
                if yroi==27
                    fisherz_roi(4,signif_overlap)=67;
                end 
                if yroi==28||yroi==29||yroi==30
                    fisherz_roi(4,signif_overlap)=7;
                end 
                
            end
        end
    end
    DFC_statistic_cal{clu_num+1,29} = fisherz_roi;
    
    
end
for clu_num = 1:nclu
    overlap_sig_r_t = DFC_statistic_cal{clu_num+1,5} + DFC_statistic_cal{clu_num+1,9};
    DFC_statistic_cal{clu_num+1,21} = overlap_sig_r_t;
end