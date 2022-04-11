%% Function calling

size = 6;
% creating a positive definite matrix
Matrix = complex(randn(size), randn(size)) ;
Matrix = Matrix'*Matrix;
% calling inbuilt inverse operation for comparing
invM = inv(Matrix);
choleskyInvM = choleskyInverse(Matrix);
fprintf("Error : %.4f\n", norm(choleskyInvM-invM));

%% Function defination
function invMat = choleskyInverse(Matrix)
    N = length(Matrix);
    Lower = zeros(N,N);
    Upper = zeros(N,N);
    Lower(1,1) = sqrt(Matrix(1,1));
    Upper(1,1) = Lower(1,1);
    for a=2:N
        Lower(a,1) = Matrix(a,1)/Lower(1,1);
        Upper(1,a) = Matrix(1,a)/Lower(1,1);
    end
    
    for i=2:N
        for j = i:N
            if i == j
                Lower(i,j) = sqrt(Matrix(j,i)-Lower(j,1:i-1)*Upper(1:i-1,i));
                Upper(j,i) = Lower(j,i);
            else
                Lower(j,i)=(Matrix(j,i)-Lower(j,1:i-1)*Upper(1:i-1,i))/Lower(i,i);
            end
        end
        for k = i+1:N
            Upper(i,k) = (Matrix(i,k)-Lower(i,1:i-1)*Upper(1:i-1,k))/Lower(i,i);
        end
    end
    
    %calculating inverse
    
    Linv = Lower \ speye(size(Lower));
    Uinv = Upper \ speye(size(Upper));
    invMat = Uinv*Linv;
end
