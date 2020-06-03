function error = costFunctionclassicK(Fs, x)

% cost function classical Kruppa 
% input FS : fundemental matrix 
% x : initial gauss for intrixsic parameter
% output: scalar which is norm of the matrix 

A=[x(1) x(2) x(3); 0 x(4) x(5); 0 0 1];

error = [];

w = A * A';

function M  = toMat( x )
    
 M=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];

end

function epi = epipole(Fs)

[U, S, V]  = svd(Fs');

epi = V(:,end);

epi = toMat(epi);

end


for i = 1 :  10 
    for j = i+1 : 10
            %left part 
            p1 = (Fs(:,:,i,j) * w * Fs(:,:,i,j)')/norm(Fs(:,:,i,j) * w * Fs(:,:,i,j)','fro');
            %right part 
            p2 = epipole(Fs(:,:,i,j)) * w * epipole(Fs(:,:,i,j))'/norm(epipole(Fs(:,:,i,j)) * w * epipole(Fs(:,:,i,j)), 'fro');
            %cost fun
            cost =norm(p1 - p2,'fro');
            
            error = [error cost];
            %error= 1000* error
        
    end
end


end

