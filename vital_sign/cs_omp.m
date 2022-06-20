%Victor Barranca, Yolanda Hu, Alex Xuan
%************************************************************************%
function hat_x=cs_omp(y,T_Mat,m)
% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

n=length(y);
s=floor(n/4);     %/4 typ    : controls sparsity and number of iterations   
%s = 750
hat_x=zeros(1,m);                                 %  empty soln                    
Aug_t=[];                                         %  set of highly correlated columns: initialize to be empty
r_n=y;                                            %  residual: initialize with rhs 

for times=1:s;                                  %  number of iterations (more, more dense)
    product=abs(T_Mat'*r_n);                 %look at correlations bn residual and cols of samp mtx
    
    [val,pos]=max(product);                       %  find maximal correlation and its corresponding column
    Aug_t=[Aug_t,T_Mat(:,pos)];                   %  augment the set of highly correlated columns
    T_Mat(:,pos)=zeros(n,1);                      %  don't consider this highly corr column in future correlation computations
    aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  solve min augt x-y
    r_n=y-Aug_t*aug_x;                            %  compute new residual
    pos_array(times)=pos;                         %  
    
end
hat_x(pos_array)=aug_x;                           %  the value for the solution estimate in component posarray_i is aug_x_i, rest entries 0, not each subsequent residual is orthog to all prev considered highly corr columns
