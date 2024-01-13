function [stokesVectors S0 LAMBDA dop] = JonesToStokes(J)
s = size(J);
lambdaCount = s(1);
N = s(2)./2;
[LAMBDA] = generateGeneralizedStokesBasis(N);
s = size(LAMBDA);
stokesCount = s(1);
stokesVectors = zeros(lambdaCount,stokesCount,'single');
N=s(2);
S0 = zeros(1,lambdaCount,'single');
for lambdaIdx=1:lambdaCount
    j = squeeze(J(lambdaIdx,:)).';
     j = j./sqrt(sum(sum(abs(j).^2)));
    S0(lambdaIdx) = sum(abs(j).^2);
    
    for stokesIdx=1:stokesCount
        % L = squeeze(LAMBDA(stokesIdx,:,:))./sqrt((N-1)./(N/2));
        L = squeeze(LAMBDA(stokesIdx,:,:));
        stokesVectors(lambdaIdx,stokesIdx) = real(sum(conj(j).*(L*j)))./sqrt(2.*(N-1)./(N));
        % size(L)
        % size(j)
        
        % for i=1:length(j)
        %    for k=1:length(j)
        %    stokesVectors(lambdaIdx,stokesIdx) = stokesVectors(lambdaIdx,stokesIdx)+conj(j(i)).*L(i,k).*j(k);
        %    end
        % end
        % stokesVectors(lambdaIdx,stokesIdx) = trace((L*(j*j'));
    end
end
stokesVectorTotal = sum(stokesVectors);
dop = sqrt(sum(abs(stokesVectorTotal).^2))./sum(S0);
end

