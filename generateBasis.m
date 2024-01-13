function [MODES,L,M,O] = generateBasis(maxMG,TH,R,mfd)
if (maxMG<0)
    lpBasis = 1;
else
    lpBasis = 0;
end
maxMG = abs(maxMG);
s = size(TH);
Nx = s(1);
Ny = s(2);

L=0:maxMG;
M=1:maxMG;
[L,M] = meshgrid(L,M);

MG = 2.*M+L-1;
MG = MG(1:end);
[mg mgIdx] = sort(MG);
idxes = find(mg<=maxMG);
L = L(mgIdx(idxes));
M = M(mgIdx(idxes));
%O = O(mgIdx(idxes));
O = L;
modeCount = length(L(L==0))+2.*length(L(L>0));

MODES = zeros(modeCount,Nx,Ny,'single');
LL = zeros(1,modeCount);
MM = zeros(1,modeCount);
OO = zeros(1,modeCount);

idx=1;
for lIdx=1:length(L)
    l = L(lIdx);
    % for mIdx=1:length(M)
    m = M(lIdx);
    
    if (lpBasis)
        [Ex,Ey] = LPmodeMFD(TH,R,mfd,l,m);
    else
        [Ex,Ey] = OAMmodeMFD(TH,R,mfd,l,m);
    end
    MODES(idx,:,:) = Ex;
    
    
    if (l==0)
        OO(idx)=-1;
        LL(idx) = l;
        MM(idx) = m;
    else
        OO(idx)=0;
                LL(idx) = l;
        MM(idx) = m;
    end
    idx=idx+1;
    
    if (l>0)
        MODES(idx,:,:) = Ey;
        OO(idx) = 1;
                LL(idx) = l;
        MM(idx) = m;
        idx=idx+1;
        
    end
  %  LL(lIdx) = l;
  %  MM(lIdx) = m;
  %  OO(lIdx) = O(idx-1);
end
L = LL;
M = MM;
O = OO;
end

