%本函数采用Jacobi方法计算实对称矩阵的所有特征值和对应特征向量 
%返回值D为特征值对角阵，V为对应特征向量构成的正交方阵 
%即有V'*A*V=D,V'*V=I 
%采用查找绝对值最大的非对角元素方法 
function [D,V]=Jaco(A) 
%     tic; 
%     %检验输入是否合法 
    b=size(A); 
%     if b(1)~=b(2)    %行列不等 
%         error('MATLAB:Jaco:Invalid Matrix,The Matrix input should be a Symmetry Phalanx.  See Jaco.'); 
%     end 
    n=max(b); 
%     for i=1:n    %非对称 
%         for j=1:n 
%             if abs(A(i,j)-A(j,i))>eps    %不能用不等号，因为机器有误差 
%                 error('MATLAB:Jaco:Invalid Phalanx,The Phalanx input should be a Symmetric one.  See Jaco.'); 
%             end 
%         end 
%     end 
     
    %实际计算 
    %初始化D，V为单位矩阵，避免多次分配空间 
    %并且相乘不会造成影响 
    D=eye(n); 
    V=eye(n); 
    %采用扫描绝对值最大A(p,q)的算法 
    p=0;    %储存最大元素所在行 
    q=0;    %储存最大元素所在列 
    maxpq=0;    %储存绝对值最大元素 
    for i=1:n-1 
        for j=i+1:n 
            if abs(A(i,j))>abs(maxpq) 
                maxpq=A(i,j); 
                p=i; 
                q=j; 
            end 
        end 
    end 
    % 
    counter = 0;
    while abs(maxpq)>eps 
        counter = counter + 1;
        maxpq=0;%务必清零，否则会死循环 
        phi=atan2(2*A(p,q),A(p,p)-A(q,q))/2; 
        U=eye(n); 
        U(p,p)=cos(phi); 
        U(q,q)=cos(phi); 
        U(p,q)=-sin(phi); 
        U(q,p)=sin(phi);
        D=U'*A*U;
        V=V*U;
        A=D;
        
        %改写maxpq
        for i=1:n-1
            for j=i+1:n
                if abs(A(i,j))>abs(maxpq)
                    maxpq=A(i,j);
                    p=i;
                    q=j;
                end
            end
        end
%         if(counter>10*n)
%             break;
%         end
    end
%     toc;
end