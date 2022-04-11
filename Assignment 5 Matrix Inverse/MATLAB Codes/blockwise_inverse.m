%% Function calling

size = 6;
Matrix = randn(size);
blocksizes =  2 * ones(1, 3);
% calling inbuilt inverse operation for comparing
invM = inv(Matrix);
blockInvM = blockInverse(Matrix, blocksizes);
fprintf("Error : %.4f\n", norm(blockInvM-invM));

%% function for block inversion

function m_inverse = blockInverse(M, framesizes, pass)

    if nargin < 3 || ~pass
        % to check whether matrix passed to function valid or not
        assert(ismatrix(M));
        assert(sum(framesizes) == size(M, 1) && size(M, 1) == size(M, 2));
    end

    if length(framesizes) == 1
        % at last recursive call
        m_inverse = M \ speye(size(M));
    else
        % getting framesize for each iterations
        frameSize = framesizes(1);            
        % Slicing matrix into four sublevel matrices
        m11 = M(1:frameSize, 1:frameSize);    
        m12 = M(1:frameSize, frameSize+1:end);
        m21 = M(frameSize+1:end, 1:frameSize);
        m22 = M(frameSize+1:end, frameSize+1:end);
        
        %applying block inversion formula
        invA = m11 \ speye(size(m11));
        Ainv_B = m11 \ m12;
        C_Ainv = m21 / m11;
        % recursive call
        n22 = blockInverse(m22 - m21 * (m11 \ m12), framesizes(2:end), true); 
        n21 = -n22 * C_Ainv;
        n12 = -Ainv_B * n22;
        n11 = invA - n12 * C_Ainv;
        
        % computed matrix inverse at each recursive calls and parent call
        m_inverse = [n11 n12; n21 n22];
    end

end
