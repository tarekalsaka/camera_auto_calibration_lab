function error = costSimpifiedK( Fs, x )
% cost function to simplified Kruppas equations
%input
%Fs:Fundamental matrices
%x:The initial instrinsic parameters
% output: coulm vector (1,45)
% if you take all the combination for substract result will get 45 *3 

A=[x(1) x(2) x(3); 0 x(4) x(5); 0 0 1];

w = A * A';

error = [];

for i = 1 : 9
    for j = i+1 : 10
        
        [U,D,V] = svd(Fs(:,:,i,j));
        
        u1 = U(:,1);
        u2 = U(:,2);
        u3 = U(:,3);
        
        v1 = V(:,1);
        v2 = V(:,2);
        v3 = V(:,3);
        
        r = D(1,1);
        s = D(2,2);
        
        %Simplified Kruppa's Equation
        A = (r^2 * v1' * w * v1) * pinv(u2' * w * u2);
        B = (r * s * v1' * w * v2) * pinv(-u1' * w * u2);
        C = (s^2 * v2' * w * v2) * pinv(u1' * w * u1);
        
        cost1  = A - B;
        cost2 = B - C;
        cost3 = C - A;
        
        error = [error cost1 ];
        % or error = [error cost1 cost2 cost3];

        
    end
end
end
