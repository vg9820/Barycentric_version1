function [M,b] = laplacian_mean(x,t,CBoundary,Circle)
%创建稀疏矩阵，用mean weight
n = size(x,1);
M = sparse(n,n);
b = zeros(n,2);
for i = 1:n
    index = find(CBoundary == i);
    %判断是否为边界点
    if(size(index)~=0)
        M(i,i) =1;
        b(i,:) = Circle(index,:);
        continue;
    end
    %内点情况,需要统计loop元素的个数
    Loop = OneRing(i,t);
    %mean weight
    omega_sum = 0;
    m = length(Loop);
    omega = zeros(m,1);
    for j = 1:m
        j1 = j;
        j0 = mod(j-1,m);
        if(j0 == 0)
            j0 = m;
        end
        j2 = mod(j+1,m);
        if(j2 == 0)
            j2 = m;
        end
        v0 = x(i,:);
        v1 = x(Loop(j1),:);
        v2 = x(Loop(j0),:);
        v3 = x(Loop(j2),:);
        e1 = v1 - v0;
        e2 = v2 - v0;
        e3 = v3 - v0;
        
        E_len  = norm(e1);
        alpha1 = acos(dot(e1,e2)/(norm(e1)*norm(e2)));
        alpha2 = acos(dot(e1,e3)/(norm(e1)*norm(e3)));
        omega(j) = (tan(0.5*alpha1)+tan(0.5*alpha2))/E_len;
        omega_sum = omega_sum + omega(j);
    end
    omega = -omega./omega_sum;
    for j = 1:m
        M(i,Loop(j)) = omega(j);
    end
    M(i,i) = 1;
end

end



