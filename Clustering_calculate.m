function Clustering_calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-means clustering
% written by
% Qi Liu; Siyu Zhu
% Weihua Zhao: zarazhao@uestc.edu.cn
% last edited Aug 2023

%input------DFC_cha_inf

%output----DFC_unfold_inf£º
%       ---Kmeans_results£ºK-means results
%       ---Kmeans_results.IDX: labels for each functional connectivity connectome
%       ---Kmeans_results.C:center for each cluster
%       ---Tol_CSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('DFC_cha_inf.mat')
nsub=82;     
nchannel=30;     
nwin=170;    %window number
maxk=10;   % max cluster number

DFC_unfold_inf = zeros(nsub*nwin,(nchannel-1)*nchannel/2+1);
for sub_num = 1:nsub
    
    for win_num = 1:nwin
        DFC_unfold_inf((sub_num-1)*nwin+win_num,1)=sub_num+200; 
        dfc_num=1;
        for  xroi = 1 : nchannel-1
            for yroi = xroi+1 : nchannel              
                DFC_unfold_inf((sub_num-1)*nwin+win_num,dfc_num+1)=DFC_cha_inf{sub_num+1,2}(xroi,yroi,win_num); 
                dfc_num=dfc_num+1;                                                                             
            end
        end
    end
    
end



for k=2:maxk  
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(DFC_unfold_inf(:,2:436),k,'Distance','cityblock','Replicates',20,'Display','final'); %,'Options',opt);   
    Kmeans_results{k}.IDX=IDX;   
    Kmeans_results{k}.C=C;       
    Kmeans_results{k}.SUMD=SUMD; 
    Kmeans_results{k}.D=D; 
end


Tol_CSS=zeros(maxk,1);    
for k=2:maxk
     disp(['Calculating for ' num2str(k) 'clusters'])
     
     for clu_num = 1:k
         for win_num = 1:13940     
             if Kmeans_results{1,k}.IDX(win_num) == clu_num
                single_CSS = (DFC_unfold_inf(win_num,2:421) - Kmeans_results{1,k}.C(clu_num,1:420))*(DFC_unfold_inf(win_num,2:421) - Kmeans_results{1,k}.C(clu_num,1:420))';
                Tol_CSS(k)=Tol_CSS(k)+single_CSS;
             end
         end
     end
end
plot(2:maxk,Tol_CSS(2:maxk));


