function error=costFunItrativeMethod()

DIAQ = IDAQlinearsolve();

global PPM; 
global A;
error = [];
for i=1 : 10
    P= PPM(:,:,i);
    
    cost= (norm((A*A'-P*DIAQ*P'),'fro'))^2   ;
    
    error= [error cost];
    
end 
end



    
    
    
    