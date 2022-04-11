%% Function calling

size = 6;
Matrix = randn(size);
% calling inbuilt inverse operation for comparing
invM = inv(Matrix);
gaussInvM = gaussInverse(Matrix);
fprintf("Error : %.4f\n", norm(gaussInvM-invM));

%% Function defination

function invMatrix =  gaussInverse(M)
    % getting size of matrix
    [rows, ~]  = size(M);
    % creating identity matrix for augumentation
    invMatrix = eye(rows);
    for j = 1 : rows
        for i = j : rows
            if M(i,j) ~= 0
                for k = 1 : rows
                    s = M(j,k); M(j,k) = M(i,k); M(i,k) = s;
                    s = invMatrix(j,k); invMatrix(j,k) = invMatrix(i,k); invMatrix(i,k) = s;
                end
                t = 1/M(j,j);
                for k = 1 : rows
                    M(j,k) = t * M(j,k);
                    invMatrix(j,k) = t * invMatrix(j,k);
                end
                for L = 1 : rows
                    if L ~= j
                        t = -M(L,j);
                        for k = 1 : rows
                            M(L,k) = M(L,k) + t * M(j,k);
                            invMatrix(L,k) = invMatrix(L,k) + t * invMatrix(j,k);
                        end
                    end
                end           
            end
            break
        end
        if M(i,j) == 0
            disp('Warning: Singular Matrix')
            invMatrix = Inf(rows);
            return
        end
    end
end