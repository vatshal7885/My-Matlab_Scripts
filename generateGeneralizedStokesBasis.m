function [LAMBDA] = generateGeneralizedStokesBasis(N)
%N = 1;
%4.*N^2-1
pauli_1 = [[1 0].',[0 -1].'].';
pauli_2 = [[0 1].',[1 0].'].';
pauli_3 = [[0 -1i].',[1i 0].'].';
pauli0 = eye(2,'single');
pauli = zeros(3,2,2,'single');
pauli(1,:,:) = pauli_1;
pauli(2,:,:) = pauli_2;
pauli(3,:,:) = pauli_3;

LAMBDA = zeros(4.*N.^2-1,2.*N,2.*N,'single');
GellMann = 1;
if (~GellMann)
    isDiagBlock = zeros(2.*N);
    dx=0;
    idx=1;
    normFctr = sqrt(N);
    for i=1:N
        for j=1:3
            LAMBDA(idx,dx+(1:2),dx+(1:2)) = squeeze(pauli(j,:,:))./normFctr;
            isDiagBlock(dx+(1:2),dx+(1:2)) = 1;
            idx=idx+1;
        end
        dx=dx+2;
    end
    
    for j=1:2.*N
        for k=1:2.*N
            if (~isDiagBlock(j,k))
                if (j<k)
                    LAMBDA(idx,j,k) = 1./normFctr;
                    LAMBDA(idx,k,j) = 1./normFctr;
                else
                    LAMBDA(idx,j,k) = -1i./normFctr;
                    LAMBDA(idx,k,j) = 1i./normFctr;
                end
                % X(j,k) = 1;
                idx=idx+1;
            end
        end
    end
    %fprintf('2 %i,%i\n',4.*N^2-N,idx-1);
    
   % normFctr = sqrt(N./((N-1).^2+(N-1)))
    for n=1:(N-1)
        dx=0;
        normFctr = sqrt(N./((n).^2+(n)))
        for i=1:n
            LAMBDA(idx,dx+(1:2),dx+(1:2)) = pauli0./normFctr;
            %idx=idx+1;
            dx=dx+2;
        end
        LAMBDA(idx,dx+(1:2),dx+(1:2)) = -n.*pauli0./normFctr;
        idx=idx+1;
    end
else
    N = 2.*N;
    idx=1;
    for k=1:(N-1)
        normFctr = sqrt(2./(k.*(k+1)));
        for j=1:k
            LAMBDA(idx,j,j) = 1.*normFctr;
        end
        LAMBDA(idx,k+1,k+1) = -k.*normFctr;
        idx=idx+1;
    end
    idx
    %for k=1:(N^2-N)
    for j=1:N
        for k=1:N
            if (k<j)
                LAMBDA(idx,j,k) = 1;
                LAMBDA(idx,k,j) = 1;
                idx=idx+1;
                LAMBDA(idx,j,k) = 1i;
                LAMBDA(idx,k,j) = -1i;
                idx=idx+1;
            end
        end
    end
    fprintf('Basis size = %i,%i\n',4.*(N/2)^2-1,idx-1);
end
