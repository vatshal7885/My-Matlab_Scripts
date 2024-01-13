includeStokes = 1;
doRandom = 0;
dop = nan(totalModeCount,lambdaShortCount);
dopLG = nan(totalModeCount,lambdaShortCount);
dopMin = zeros(1,totalModeCount,'single');
dopLGMin = zeros(1,totalModeCount,'single');
[minV lambda0Idx] = min(abs(lambdaShort-1545));
bandwidth = nan(lambdaShortCount,1);
c = 299792458;
thresh = 10.^(-1./10)

for modeIdx=1:totalModeCount
    for dLambda=1:(lambdaShortCount)
        idx1 = lambda0Idx-dLambda;
        idx2 = lambda0Idx+dLambda;
        if (idx1>=1 & idx2<=lambdaShortCount)
            dop(modeIdx,dLambda) = sum(sum(abs(ResponseSpectral_PS(lambda0Idx+(-dLambda:dLambda),modeIdx)).^2))./length(-dLambda:dLambda);
            dopLG(modeIdx,dLambda) = sum(sum(abs(ResponseSpectral_LG(lambda0Idx+(-dLambda:dLambda),modeIdx)).^2))./length(-dLambda:dLambda);
            f1 = c./(lambdaShort(lambda0Idx-dLambda).*1e-9);
            f2 = c./(lambdaShort(lambda0Idx+dLambda).*1e-9);
            bandwidth(dLambda) = 1e-12.*abs(f1-f2);
        end
    end
    dopLGMin(modeIdx) = sum(sum(abs(ResponseSpectral_LG(:,modeIdx)).^2))./lambdaShortCount;
    dopMin(modeIdx) = sum(sum(abs(ResponseSpectral_PS(:,modeIdx)).^2))./lambdaShortCount;
end


cutoff = ones(1,totalModeCount).*max(bandwidth);
cutoffLG = ones(1,totalModeCount).*max(bandwidth);
%f = c./(lambdaShort.*1e-9);
for modeIdx=1:totalModeCount
    foundIt = 1;
    foundItLG = 1;
    for dLambda=1:lambdaShortCount
        if (dop(modeIdx,dLambda)<thresh & foundIt)
            cutoff(modeIdx) = bandwidth(dLambda);
            foundIt = 0;
        end
        if (dopLG(modeIdx,dLambda)<thresh & foundItLG)
            cutoffLG(modeIdx) = bandwidth(dLambda);
            foundItLG = 0;
        end
    end
end

%Random input
if (doRandom)
[minV lambda0IdxR] = min(abs(lambda-1545));
randCount = 10;
ResponseSpectral_Random = zeros(lambdaShortCount,randCount,'single');
for randIdx=1:randCount
    v = rand(1,totalModeCount,'single').*exp(1i.*2.*pi.*rand(1,totalModeCount,'single'));
    v = v./sqrt(sum(sum(abs(v).^2)));
    v = v.';
    
    Mc = squeeze(MTMs(lambda0IdxR,:,:));
    v0 = conj((Mc*v));
    v0 = v0./sqrt(sum(sum(abs(v0).^2)));
    for lambdaIdx=1:lambdaShortCount
        M0 = squeeze(MTMs((lambdaIdx-1).*lambdaStep+1,:,:));
        v1 = M0*v;
        v1 = v1./sqrt(sum(sum(abs(v1).^2)));
        ResponseSpectral_Random(lambdaIdx,randIdx) = abs(sum(v0.*(v1))).^2;
    end
end

dopR = nan(randCount,lambdaShortCount);
for modeIdx=1:randCount
    for dLambda=1:(lambdaShortCount)
        idx1 = lambda0Idx-dLambda;
        idx2 = lambda0Idx+dLambda;
        if (idx1>=1 & idx2<=lambdaShortCount)
            dopR(modeIdx,dLambda) = sum(sum(abs(ResponseSpectral_Random(lambda0Idx+(-dLambda:dLambda),modeIdx)).^2))./length(-dLambda:dLambda);
        end
    end
end

cutoffR = ones(1,randCount).*max(bandwidth);
%f = c./(lambdaShort.*1e-9);
for modeIdx=1:randCount
    foundIt = 1;
    for dLambda=1:lambdaShortCount
        if (dopR(modeIdx,dLambda)<thresh & foundIt)
            cutoffR(modeIdx) = bandwidth(dLambda);
            foundIt = 0;
        end
    end
end

figure(2342341);
subplot(1,2,1);
plot(ResponseSpectral_Random);
subplot(1,2,2);
plot(bandwidth,dopR.');
xlim([0 5]);
end
%End Random input


figure(324234234);
subplot(1,3,1);
plot(bandwidth,dop.');
ylim([0 1]);
grid on;
ylim([0 1]);
xlim([0 5]);
xlabel('\Deltaf (THz)');
axis square;

subplot(1,3,2);
plot(bandwidth,dopLG.');
ylim([0 1]);
xlim([0 5]);
xlabel('\deltaf (THz)');
grid on;
axis square;

subplot(1,3,3);
x0 = zeros(1,8);for i=1:8 x0(i) = 2.*sum(1:i);end
hold off;
plot(1:totalModeCount,cutoffLG,'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none','Color',[1 0 0]);
hold on;
plot(1:totalModeCount,cutoff,'MarkerFaceColor',[0 0.5 0],'Marker','o','LineStyle','none','Color',[0 0.5 0]);
plot([0 totalModeCount+1],max(bandwidth).*[1 1],'-black');
%axis square;
ylim([0 max(bandwidth).*1.1]);
xlim([0 totalModeCount+1]);
set(gca,'YTick',0:6);
set(gca,'XTick',x0+0.5);
grid on;

figure(2342342);
plot(1:totalModeCount,dopMin,'x',1:totalModeCount,dopLGMin,'o');
ylim([0 1]);
%pause;
purity = zeros(1,totalModeCount,'single');
if (includeStokes)
    dopStokesLG = zeros(1,totalModeCount);
    dopStokes = zeros(1,totalModeCount);
    for modeIdx=1:totalModeCount
        D = 0;
        Response = zeros(lambdaCount,totalModeCount,'single');
        ResponseLG = zeros(lambdaCount,totalModeCount,'single');
        Pin = squeeze(PSPin(:,modeIdx));
        for lambdaIdx=1:lambdaCount
            J = squeeze(MTMs(lambdaIdx,:,:));
            Pout = J*Pin;
            Pout = Pout./sqrt(sum(sum(abs(Pout).^2)));
            d = ((Pout)*(Pout'))./lambdaCount;
            D = D+d;
            Response(lambdaIdx,:) = Pout;
            ResponseLG(lambdaIdx,:) = squeeze(J(:,modeIdx));
        end
        p = zeros(totalModeCount,totalModeCount,'single');
        for m=1:totalModeCount
            u_m = zeros(1,totalModeCount,'single');
            u_m(m) = 1;
            for n=1:totalModeCount
                u_n = zeros(1,totalModeCount,'single');
                u_n(n) = 1;
                p(m,n) = u_m*D*(u_n');
            end
        end
        purity(modeIdx) = trace(p*p);
        fprintf('%i %3.3f\n',modeIdx,purity(modeIdx));
        [stokesVectors S0 LAMBDA dops] = JonesToStokes(Response);
        dopStokes(modeIdx) = dops;        
    end
    
    figure(1);
    plot(1:totalModeCount,real(purity),'o',1:totalModeCount,dopStokes,'x',1:totalModeCount,dopMin,'+');
    
    figure(2142342);
    subplot(1,3,3);

    hold off;
    plot(1:totalModeCount,squeeze(spatialCoherenceLGStates(2,1,:)),'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none','Color',[1 0 0]);
    hold on;
    plot(1:totalModeCount,squeeze(spatialCoherenceStates(2,1,:)),'MarkerFaceColor',[0 0.5 0],'Marker','o','LineStyle','none','Color',[0 0.5 0]);
    plot(1:totalModeCount,squeeze(dopStokes),'MarkerFaceColor',[0 0.5 0],'Marker','x','LineStyle','none','Color',[0 0.5 0]);
    plot(1:totalModeCount,squeeze(dopStokesLG),'MarkerFaceColor',[1 0 0],'Marker','x','LineStyle','none','Color',[1 0 0]);
    hold off;
    ylim([0 1]);
    xlim([0 totalModeCount+1]);
    axis square;
    
end

figure(330299);
x0 = zeros(1,8);for i=1:8 x0(i) = 2.*sum(1:i);end
hold off;
semilogy(1:totalModeCount,cutoffLG,'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none','Color',[1 0 0]);
hold on;
semilogy(1:totalModeCount,cutoff,'MarkerFaceColor',[0 0.5 0],'Marker','o','LineStyle','none','Color',[0 0.5 0]);
semilogy([0 totalModeCount+1],mean(cutoffR).*[1 1],'MarkerFaceColor',[0 0 1],'LineStyle','-','Color',[0 0 1]);
plot([0 totalModeCount+1],max(bandwidth).*[1 1],'-black');

xlim([0 totalModeCount+1]);

set(gca,'XTick',x0+0.5);
grid on;

