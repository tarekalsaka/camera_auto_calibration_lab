function finalQ= IDAQlinearsolve()
% function to compute the Dual Image Absolute Quadtic matrix 
% output Q 4 *4 Dual Image Absolute Quadtic matrix 

global PPM;

% Symbolic variables

syms X11 X12 X13 X14 X21 X22 X23 X24 X31 X32 X33 X34 X41 X42 X43 X44;

% Q = [X11, X12, X13, X14;
%     X21, X22, X23, X24;
%     X31, X32, X33, X34;
%     X41, X42, X43, X44];

% symitric 
Q = [X11, X12, X13, X14;
    X12, X22, X23, X24;
    X13, X23, X33, X34;
    X14, X24, X34, X44];

% Q2 = [X11, X12, X13, X14 X21, X22, X23, X24 X31, X32, X33, X34 X41, X42, X43, X44];
Q2 = [X11 X12 X13 X14 X22 X23 X24 X33 X34 X44];

for i=1:10
    P = PPM(:, :, i);

    % eqution from constraint x0, y0 = 0 
    DIAC(3*(i-1)+1:3*i,:) = P * Q * P';
    L(2*(i-1)+1:2*i,1) = DIAC(3*(i-1)+1:3*i-1,3)==0;

end
 
% size(L)
[mat,b] = equationsToMatrix(L);
% eq=mat*Q2'==0 ;
% S = solve(eq)
% %linsolve(mat,b)
[~,D,V]= svd(mat);
Lastcolv= V(:,end);
finalQ =[Lastcolv(1), Lastcolv(2), Lastcolv(3),Lastcolv(4);
         Lastcolv(2), Lastcolv(5), Lastcolv(6),Lastcolv(7);
         Lastcolv(3), Lastcolv(6), Lastcolv(8), Lastcolv(9);
         Lastcolv(4), Lastcolv(7), Lastcolv(9), Lastcolv(10)];
finalQ= double(finalQ); 

% enforc rank 3 to the matrix final Q by assgin zero to the smallest eigenvalue
[U,D,V] = svd(finalQ);
D(4,4)=0;
%retraive the Q 
finalQ= U*D*V;



