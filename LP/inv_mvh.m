% function Ad_eye=inv_mvh(B_inv,Ad)
% Ad_eye=B_inv\Ad;% ����
function A_inv = inv_mvh(A)
% �Ծ�����г����б任������

[row, col] = size(A);

% BΪ��λ����
B = eye(row);

for i = 1 : row
    % ���ν��Խ��е�Ԫ�ع�һ��
    div_i = A(i, i);
    for j = 1 : col
       A(i, j) = A(i, j) / div_i;
       B(i, j) = B(i, j) / div_i;
    end

    for ii = 1 : row
       tmp_ii = - A(ii, i) / A(i, i);
       if i == ii
           tmp_ii = 0;
       end
       
       % �����б任
       for jj = 1 : col
           A(ii, jj) = A(ii, jj) + tmp_ii * A(i, jj);
           B(ii, jj) = B(ii, jj) + tmp_ii * B(i, jj);
       end
       
    end
end

A_inv = B;