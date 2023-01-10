%����������Jacobi��������ʵ�Գƾ������������ֵ�Ͷ�Ӧ�������� 
%����ֵDΪ����ֵ�Խ���VΪ��Ӧ�����������ɵ��������� 
%����V'*A*V=D,V'*V=I 
%���ò��Ҿ���ֵ���ķǶԽ�Ԫ�ط��� 
function [D,V]=Jaco(A) 
%     tic; 
%     %���������Ƿ�Ϸ� 
    b=size(A); 
%     if b(1)~=b(2)    %���в��� 
%         error('MATLAB:Jaco:Invalid Matrix,The Matrix input should be a Symmetry Phalanx.  See Jaco.'); 
%     end 
    n=max(b); 
%     for i=1:n    %�ǶԳ� 
%         for j=1:n 
%             if abs(A(i,j)-A(j,i))>eps    %�����ò��Ⱥţ���Ϊ��������� 
%                 error('MATLAB:Jaco:Invalid Phalanx,The Phalanx input should be a Symmetric one.  See Jaco.'); 
%             end 
%         end 
%     end 
     
    %ʵ�ʼ��� 
    %��ʼ��D��VΪ��λ���󣬱����η���ռ� 
    %������˲������Ӱ�� 
    D=eye(n); 
    V=eye(n); 
    %����ɨ�����ֵ���A(p,q)���㷨 
    p=0;    %�������Ԫ�������� 
    q=0;    %�������Ԫ�������� 
    maxpq=0;    %�������ֵ���Ԫ�� 
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
        maxpq=0;%������㣬�������ѭ�� 
        phi=atan2(2*A(p,q),A(p,p)-A(q,q))/2; 
        U=eye(n); 
        U(p,p)=cos(phi); 
        U(q,q)=cos(phi); 
        U(p,q)=-sin(phi); 
        U(q,p)=sin(phi);
        D=U'*A*U;
        V=V*U;
        A=D;
        
        %��дmaxpq
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