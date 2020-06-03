clc 
clear 
load 'data.mat'
global PPM
%find the DIAQ 
Q = IDAQlinearsolve();
rank(Q)
convertQToDouble= double(Q);
det(convertQToDouble)

%Checkif Q respect the proprties det Q =0 
if det(convertQToDouble)==0
    disp("det = 0 and the rank is ")
else
    disp("det is not zero and rank is :")
    disp(rank(convertQToDouble))
    % inforce for det = 0 
    [U,D,V] = svd(Q);
    D(4,4)=0;
    %retraive the Q 
    Q2= U*D*V;
end 
d= eig(Q2);
isposdef = all(d) > 0;
%check if the Q2 positive semi-definite 
if isposdef == all(d) > 0
    disp("after endforc the det = 0 tha rankd become")
    disp(rank(Q2))
    disp(" and Q is positive semi-definite")
else 
    disp("Q is not positive semi-definite")
    
end 
% find the plane at infinty from the null space of Q 
[~,~,v]=svd(Q2);
planeAtInf= v(1:3,end)

% compute w 
P = PPM(:,:,9);

w = P * Q2 * P'

%check if w is simi positive by try to find Cholesky factorization 
try chol(w)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end
