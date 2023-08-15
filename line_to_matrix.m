% ������ת��Ϊ����
% matrix = line_to_matrix(data,nroi,diagonal) 
% input -- data  ������
%       -- nroi  roi����  ����߳�
%       -- diagonal  �Խ�����Ϊ�Ǹ���
% output -- mrtrix  �Խ���Ϊdiagonal �߳�Ϊnroi�ľ���
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