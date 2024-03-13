% function Ad_eye=inv_mch(B_inv,Ad)
% Ad_eye=B_inv\Ad;% 化简
function A_inv = inv_mch(A,row, col)
% 对矩阵进行初等行变换求其逆

% [row, col] = size(A);

% B为单位矩阵
B = eye(row);

for i = 1 : row
    % 依次将对角行的元素归一化
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
       
       % 初等行变换
       for jj = 1 : col
           A(ii, jj) = A(ii, jj) + tmp_ii * A(i, jj);
           B(ii, jj) = B(ii, jj) + tmp_ii * B(i, jj);
       end
       
    end
end

A_inv = B;