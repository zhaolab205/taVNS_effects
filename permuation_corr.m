function per_p = permuation_corr(data,label,n_per,obs_r)
per_corr = zeros(n_per,1);
for  per_num = 1 : n_per
    per_index = randperm(length(label));
    per_label = label(per_index);
    [per_corr(per_num),p]=corr(data,per_label);
end
length_p = min(length(find(per_corr>obs_r))+1,length(find(per_corr<obs_r))+1);
per_p = (length_p)/(n_per+1);