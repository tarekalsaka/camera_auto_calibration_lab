function error = CostFunctionMC(Fs,x)

% cost function to Mendonca & Cipolla 
% input FS : fundemental matrix 
% x : initial gauss for intrixsic parameter

error=[];

A=[x(1) x(2) x(3); 0 x(4) x(5); 0 0 1];

for i=1 :9
    for j =i+1 : 10 
        E= A' * Fs(:,:,i,j)* A ;
        [U,S,V] = svd(E);
        sigma1 = S(1,1);
        sigma2 = S(2,2);
        cost= 0.0222222222222222*(sigma1 - sigma2)/sigma2;
        error=[error cost];
    end 
end
end 

