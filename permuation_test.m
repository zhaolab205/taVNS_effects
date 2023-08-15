function per_p = permuation_test(data,size_one,size_two,n_per,obs_p)
per_ttest = zeros(n_per,1);
per_ttest_num = 0;
for  per_num = 1 : n_per
    per_index = randperm(size_one+size_two);
    per_data = data(per_index);
    [t,per_ttest(per_num)]=ttest2(per_data(1:size_one),per_data(size_one+1:size_one+size_two));
    if per_ttest(per_num) < obs_p
        per_ttest_num = per_ttest_num +1;
    end
    
end
per_p = (per_ttest_num+1)/(n_per+1);