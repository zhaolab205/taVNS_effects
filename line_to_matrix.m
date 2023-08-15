% 行向量转换为矩阵
% matrix = line_to_matrix(data,nroi,diagonal) 
% input -- data  行向量
%       -- nroi  roi数量  矩阵边长
%       -- diagonal  对角线置为那个数
% output -- mrtrix  对角先为diagonal 边长为nroi的矩阵
function matrix = line_to_matrix(data,nroi,diagonal) 
data_num=0;
for x_roi = 1 :nroi
    matrix(x_roi,x_roi)=diagonal;
    for y_roi = x_roi+1:nroi
        data_num = data_num+1;
        matrix(x_roi,y_roi)=data(data_num);
        matrix(y_roi,x_roi)=data(data_num);
    end
end