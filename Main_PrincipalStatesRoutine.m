clear all;
%close all;
%The mode-group/s we're interested in looking at
modeGroupOfInterest = [2];

%Do Degree of Coherence/Generalized Stokes analysis?
%Could be too computationally intensive for high mode-counts.
doDOP = 1;

%The set of wavelengths to calculate the principal states for
lambdaCentres = [1526 1545];
lambdaCentresCount = length(lambdaCentres);

%The wavelength to use as the reference centre wavelength for comparison of
%second-order effects.
lambdaCentreWavelength = [1545];
    
%Re adds the chromatic dispersion that the SWI removes automatically.
addChromaticDispersion = 0;

%Temporally filters the data, by only including the response within a
%window of +/- temporalWindowSize (ps) around the centre of mass for each
%mode-group. Each mode-group has it's own time window. Setting to 0 or very
%large will turn this temporal filtering off.
temporalWindowSize =20;

%Wavelength step-size in nm used to calculate principal states
%(lambda_0+/-dLAMBDA)
dLAMBDA =1;
%These principal states will be ignored for many graphing purposes (e.g.
%the fundamental modes);
ignoreModes = [ ];

%Do you want to see what the principals states look like?
visualiseStates = 0;

%Removes a common phase reference from input principal
%states calculated at different wavelengths. Not necessary, unless you
%wanted to average principal states calculated at different wavelengths.
removeInputPrincipalStatePhaseDrift = 0;

%Do you want to output the graphs to images/video?
outputToVideo = 0;
if (outputToVideo)
  videoOutputFolder = 'D:\Frames\';
  videoRes = 300;  
end

fprintf('Loading transfer matrices...');
load('Fibre01_6MF_2014_10_16.mat');
time = linspace(0,6.11,length(lambda)).*1e-9./(32768./length(lambda));
time = time.';
fprintf('Transfer matrices loaded.');

%speed of light, you should know that one.
c0 = 299792458;
%Optical frequency
w = 2.*pi.*c0./(lambda.*1e-9);

%Total number of wavelength points
lambdaCount = length(lambda);

%Includes coupling between mode-groups. Generally better to not include
%such coupling as this is introduced by the SLM mux/demux rather than
%genuine coupling in the fibre itself. Effectively this spatially filters
%the measurement data.
includeIntraGroupCoupling = 0;

%Measure of the average spatial coherence
spatialCoherence = nan(lambdaCentresCount,length(dLAMBDA));
%Spatial coherence of the worst mode
spatialCoherenceWorst = nan(lambdaCentresCount,length(dLAMBDA));

%Forces the data to take a unitary form by performing the singular value
%decomposition, then reconstituting the matrix with U*V' (eigenvalues set
%to 1). Can have a filtering effect on the data, but could also be useful
%for theoretical work where realistic values are required, than nonetheless
%exhibit ideal behaviour.
forceUnitary = 0;

%Fits the phase of the first mode in, first mode out to a parabola, and
%takes that phase off the other modes. Removes chromatic dispersion common
%to all modes. Could be dicey if the first mode isn't the fundamental.
useFirstModeAsChromaticDispersionReference = 0;

%Decreases the wavelength resolution by a factor of lambdaStep.
%Increasing lambdaStep decreases the computation time.
lambdaStep = 8;
lambdaIdxes = 1:lambdaStep:lambdaCount;
lambdaShort = lambda(lambdaIdxes);
lambdaShortCount = length(lambdaShort);
fShort = c0./(lambdaShort.*1e-9);
N = length(lambdaShort);
dt = N/(max(fShort)-min(fShort));
timeShort =dt.*([-N/2:N/2-1]/N);
fprintf('Max Delays given lambdaStep = +/-%3.3f (ps)\n',max(timeShort).*1e12);

if (addChromaticDispersion)
    L_SMF28 = 10+2+2+2;
    L_DUT = 116.1659-L_SMF28;
    CD_SMF28 = -18;
    CD = CD_SMF28.*(L_DUT-L_SMF28)./1000 %ps/nm
    
    GV_ps = CD.*linspace(0,1,lambdaCount).*(lambda(end)-lambda(1));
    
    GV = GV_ps.*1e-12;
    CD_PHI = zeros(1,lambdaCount);
    for lambdaIdx=1:(lambdaCount-1)
        dw = w(lambdaIdx+1)-w(lambdaIdx);
        CD_PHI(lambdaIdx+1) = CD_PHI(lambdaIdx)+GV(lambdaIdx).*dw;
    end
    CD_PHI = exp(1i.*CD_PHI).';
end

%GRAPHING OPTIONS
%Figure offset for enumerating Supplementary figures;
SupFigOffset = 1000;
%Tick spacing on the dB axis
yTicks = -21:3:0;
%Tick spacing on the wavelength axis.
lambdaTicks = 1525:5:1565;
%Text labels for the wavelength axis
lambdaLabels = cell(1,length(lambdaTicks));
%Only show the wavelength for every second tick
for i=1:2:length(lambdaLabels)
    lambdaLabels(i) = cellstr(num2str(lambdaTicks(i)));
    lambdaLabels(i+1) = cellstr('');
end
%Width of lines in graph
lineWidth = 1;
%Time unit, ps
TUNITS = 1e12;
%Wavelength axis (nm)
XLIM = 1545+20.*[-1 1];
%dB axis
YLIM = [-20 1];
%Time axis (ps)
TLIM = 4.*[-1 1];

%Pixel resolution of visualised states
Npixels = 64;
%Mode-field diameter of the visualised states (normalised to the width of
%the image).
mfd = 1;
%Colormap to visualise the states with
CMP_JET = jet(256);
%Grayscale colormap used for plotting
CMP_GRAY = gray(256);
%Rotate the visualised states by dTH amount (to line up with the CCD camera
%which is at an angle (11 deg) relative to the polarization axis of the SLM)
dTH = 11.*pi./180;

%END GRAPHING OPTIONS

%Pads out the L,M indices of the modes to include polarization.
l = zeros(1,2.*length(L));
m = l;
l(1:2:end) = L;
l(2:2:end) = L;
m(1:2:end) = M;
m(2:2:end) = M;
%Calculate near-degenerate mode-group number
MG = 2.*m+l-1;
MGtick = zeros(1,max(MG)-1);
for i=1:(length(MG)-1)
    if (MG(i)~=MG(i+1))
        MGtick(MG(i)) = i+0.5;
    end
end

%MG = [1 1 2 2 2 2 3 3 4 4 4 4];
if (useFirstModeAsChromaticDispersionReference)
    unwrapIt = unwrap(angle(squeeze(MTMs(:,1,1))));
    P = polyfit(lambda(1:16:end),unwrapIt(1:16:end),2);
    p = polyval(P,lambda);
end

%Pick off the mode-group you're interested in
MGidxes = nan(1,length(MG),'single');
idx = 0;
for i=1:length(modeGroupOfInterest)
    for j=1:length(MG)
        if (MG(j)==modeGroupOfInterest(i))
            idx=idx+1;
            MGidxes(idx)=j;
            
        end
    end
end
s = size(MTMs);
%The total number of spatial/polarization modes
totalModeCount = s(2);
if (idx>s(2))
    idx = totalModeCount;
end
MGidxes = MGidxes(1:idx);
MTMs = MTMs(:,MGidxes,MGidxes);
MG = MG(MGidxes);

s = size(MTMs);
%The total number of spatial/polarization modes
totalModeCount = s(2);
spatialCoherenceStates = zeros(lambdaCentresCount,length(dLAMBDA),totalModeCount,'single');
spatialCoherenceLGStates = zeros(lambdaCentresCount,length(dLAMBDA),totalModeCount,'single');
%Colors to be used for plotting
lineColor = lines(totalModeCount);

if (outputToVideo)
    videoObjsIn = cell(1,totalModeCount);
    videoObjsOut = cell(1,totalModeCount);
    fps = 25;
    for i=1:totalModeCount
        fName = sprintf('%sPSPin_%2.2i.avi',videoOutputFolder,i);
        writerObjResponse = VideoWriter(fName,'Uncompressed AVI');
        writerObjResponse.FrameRate = fps;
        open(writerObjResponse);
        videoObjsIn(i) = {writerObjResponse};
        
        fName = sprintf('%sPSPout_%2.2i.avi',videoOutputFolder,i);
        writerObjResponse = VideoWriter(fName,'Uncompressed AVI');
        writerObjResponse.FrameRate = fps;
        open(writerObjResponse);
        videoObjsOut(i) = {writerObjResponse};
        
    end
    
end

if (includeIntraGroupCoupling==0)
    for i=1:totalModeCount
        for j=1:totalModeCount
            if (MG(i)~=MG(j))
                MTMs(:,j,i) = 0;
            end
        end
    end
end
%Puts chromatic dispersion back in, that was removed by the SWI. Not
%necessary unless you care about absolute chromatic dispersion.
if (addChromaticDispersion)
    for i=1:totalModeCount
        for j=1:totalModeCount
            MTMs(:,j,i) = MTMs(:,j,i).*CD_PHI;
        end
    end
end

%Removes chromatic dispersion common to all modes, assuming the first mode
%in/out can be used as a reference.
if (useFirstModeAsChromaticDispersionReference)
    for i=1:totalModeCount
        for j=1:totalModeCount
            MTMs(:,j,i) = MTMs(:,j,i).*exp(-1i.*p);
        end
    end
end

%Builds images of the mode basis
if (visualiseStates || outputToVideo)
    maxMG = max(MG);
    x = single(linspace(-1,1,Npixels));
    y = x;
    [X Y] = meshgrid(x,y);
    [TH R] = cart2pol(X,Y);

    TH = TH+dTH;

    [MODES L M O] = generateBasis(maxMG,TH,R,mfd);
    MGidxx = MGidxes(2:2:end)/2;
    MODES = MODES(MGidxx,:,:);
    L = L(MGidxx);
    M = M(MGidxx);
    O = O(MGidxx);
end

if (temporalWindowSize)
    window = zeros(length(time),1,'single');
    %Calculate the average total impulse response and the average impulse
    %response of a particular mode-group
    impulseTotal = 0;
    dt = 0;
    hanningWindow = ((0.5.*(1-cos((2.*pi.*(0:(lambdaCount-1)))./(lambdaCount-1)))).');
    maxMG = max(MG);
    impulseTotalMG = zeros(max(MG),length(time),'single');
    for modeIdx=1:totalModeCount
        for modeIdy=1:totalModeCount
            mg = MG(modeIdy);
            J0 = squeeze(MTMs(:,modeIdy,modeIdx));
            J0 = J0.*hanningWindow;
            t0 = abs(fftshift(fft(fftshift(J0)))./sqrt(lambdaCount)).^2;
            impulseTotal = impulseTotal+t0;
            impulseTotalMG(mg,:) = impulseTotalMG(mg,:)+t0.';
        end
    end
    %Normalize to a max of 1.
    impulseTotal = impulseTotal./max(max(impulseTotal));
    impulsePositionMG = zeros(1,max(MG));
    %Take the centre of mass of the impulse response of each mode-group
    for mgIdx=1:maxMG
        iMG = squeeze(impulseTotalMG(mgIdx,:)).';
        com = sum(iMG.*time)./sum(iMG);
        impulsePositionMG(mgIdx) = com;
    end
    tMIN = 1.1*10^12;
    tMAX = -tMIN;
    %Filter each mode-group in it's own window of size windowSize ps.
    for modeIdx=1:totalModeCount
        for modeIdy=1:totalModeCount
            mg = MG(modeIdx);
            J0 = squeeze(MTMs(:,modeIdy,modeIdx));
            j0 = fftshift(fft(fftshift(J0)));
            window = zeros(size(time),'single');
            window(abs((time-impulsePositionMG(mg)).*1e12)<temporalWindowSize) = 1;
            windowLeft = TUNITS.*(impulsePositionMG(mg))+temporalWindowSize;
            windowRight = TUNITS.*(impulsePositionMG(mg))-temporalWindowSize;
            if (windowLeft>tMAX)
                tMAX = windowLeft;
            end
            if (windowRight<tMIN)
                tMIN = windowRight;
            end
            %window = zeros(size(time),'single');
           % window(time>-175e-12) =0 ;
            j0 = j0.*window;
            J0 = ifftshift(ifft(ifftshift(j0)));
            MTMs(:,modeIdy,modeIdx) = J0;
        end
    end
    %window(time<100) = 0;
end

%View the Mode Dependent Loss
[MDL sigmas] = calcMDL(MTMs);
meanMDL = mean(10.^(MDL(20:(end-20))./10));
meanMDL = 10.*log10(meanMDL);

figure(1);
subplot(1,2,1);
plot(lambda,sigmas);
axis square;
xlabel('\lambda (nm)');
ylabel('Singular Values');
xlim(XLIM);
subplot(1,2,2);
plot(lambda,MDL);
ylim(meanMDL+3.*[-1 1]);
axis square;
xlabel('\lambda(nm)');
ylabel('MDL (dB)');
xlim(XLIM);
title(sprintf('MDL (dB) = %3.3f',meanMDL));
grid on;
%Forces the data to take a unitary form. 
if (forceUnitary)
    for lambdaIdx=1:lambdaCount
        J = squeeze(MTMs(lambdaIdx,:,:));
        [U S V] = svd(J);
        J = U*V';
        MTMs(lambdaIdx,:,:) = J;
    end
end

for dLambdaIdx = 1:length(dLAMBDA)
    dLambda = dLAMBDA(dLambdaIdx);

    [minV lambdaCentresIdx0] = min(abs(lambdaCentres-lambdaCentreWavelength));

    %array for storing the calculated in/out principal states for each
    %centre wavelength.
    PSPsout = zeros(lambdaCentresCount,totalModeCount,totalModeCount,'single');
    PSPsin = zeros(lambdaCentresCount,totalModeCount,totalModeCount,'single');
    
    GroupDelay_Eigenvalues = nan(lambdaCentresCount,totalModeCount);
    GroupDelay_CentreOfMass = nan(lambdaCentresCount,totalModeCount);
    GroupDelay_Eigenvalues_dw = zeros(lambdaCentresCount,1);
    
    for lambdaCentresIdx=1:(lambdaCentresCount);
        lambda0 = lambdaCentres(lambdaCentresIdx);

        
        %Find the indices of the wavelengths closest to that specified by
        %the user.
        [minV lambda0Idx] = min(abs(lambda0-lambda));
        [minV lambdaPlusIdx] = min(abs((lambda0+dLambda)-lambda));
        [minV lambdaMinusIdx] = min(abs((lambda0-dLambda)-lambda));

        %The three matrices that will be used to calculate the principal
        %states. Jc, centre wavelength and +/- dLambda either side.
            J0 = squeeze(MTMs(lambdaMinusIdx,:,:));
            J1 = squeeze(MTMs(lambdaPlusIdx,:,:));
            Jc  = squeeze(MTMs(lambda0Idx,:,:));
            
            %Find the maximum amplitude of all three matrices and use that
            %values as the maximum for the colorbar off all three graphs
            %below.
            maxV = max([max(max(abs(Jc))),max(max(abs(J0))),max(max(abs(J1)))]);
            
            figure(3);
            subplot(1,3,1);
            imagesc(abs(J0));
            axis square;
            axis off;
            set(gca,'CLim',[0 maxV]);
            title('J_-');
            
            subplot(1,3,2);
            imagesc(abs(Jc));
            axis square;
            axis off;
            set(gca,'CLim',[0 maxV]);
            title('J_0');
            
            subplot(1,3,3);
            imagesc(abs(J1));
            axis square;
            axis off;
            set(gca,'CLim',[0 maxV]);
            title('J_+');
            colormap(CMP_GRAY);

        %Calculate the principal states
        J0inv = J0^-1;
        T = double(J1.*(J0inv));
        [Ve De] = eigs(T,length(T),1);
        %Group delays
        GD = diag(De);
        %Output principal states
        PSPout = Ve;
        %Back-propagate through Jc to get the input principal states
        PSPin = (Jc^-1)*PSPout;

        %Keep track of the principal states from one centre wavelength to
        %the next. The solutions provided by eigs above will contain phase wraps etc, but mean they won't be consistent from 1 wavelength to the next.
        %It does this by comparing the current principal state to
        %principal states over the past 2*dLambda.
        %Indices used to track and sort the principal states
        idxes = 1:totalModeCount;
        %If this is not the first wavelength
        if (lambdaCentresIdx>1)
            %Matrix which stores the correlation between the current
            %principal states and previous wavelength principal states.
            %Uses this to try to work out which principal states at this
            %wavelength correspond to which principal states at previous
            %wavelengths.
            C = zeros(totalModeCount,totalModeCount);
            for j=1:totalModeCount
                dd = 0;
                bestIdx = 1;
                %How many centre wavelengths have the principal states been
                %calculated for in the past 2*dLambda nm?
                [minV minIdx] = min(abs(2*dLambda-abs(lambdaCentres(lambdaCentresIdx)-lambdaCentres(1:(lambdaCentresIdx-1)))));
                for k=minIdx:(lambdaCentresIdx-1)
                    p2 = (conj(squeeze(PSPsout(k,:,j)))).';
                    for i=1:totalModeCount
                        p1 = squeeze(PSPout(:,i));
                        %correlation between i th principal state and j th
                        %principal state from centre wavelength k
                        d = abs(sum(p1.*p2)).^2;
                        %Average the correlations over the 2*dLambda
                        %wavelength band.
                        C(j,i) = C(j,i)+d;
                    end
                end
            end
            
            %Use the best match for each principal state.
            for J=1:totalModeCount
                maxV = 0;
                maxIdx = 1;
                maxJdx = 1;
                
                for i=1:totalModeCount
                    for j=1:totalModeCount
                        if (abs(C(j,i))>maxV)
                            maxV = abs(C(j,i));
                            maxIdx = i;
                            maxJdx = j;
                        end
                    end
                end
                
                C(maxJdx,:) = 0;
                C(:,maxIdx) = 0;
                idxes(maxJdx) = maxIdx;
            end
        %else if this is the first centre wavelength, then just sort by
        %group-delay eigenvalues.
        else
            if (lambdaCentresIdx==1)
                [Y sortIdx] = sort(angle(GD));
                GD = GD(sortIdx);
                PSPout = PSPout(:,sortIdx);
                PSPin = PSPin(:,sortIdx);
            end
        end
        
        P_in = PSPin;
        P_out = PSPout;
        
        PSPout = PSPout(:,idxes);
        PSPin = PSPin(:,idxes);        
        GD = GD(idxes);
        
        PSPsout(lambdaCentresIdx,:,:) = PSPout;
        PSPsin(lambdaCentresIdx,:,:) = PSPin;
        
        %Locks input principal states at all wavelengths to the same phase
        if (removeInputPrincipalStatePhaseDrift)
            if (lambdaCentresIdx>1)
                for j=1:totalModeCount
                    PSPsin(lambdaCentresIdx,:,j) = PSPsin(lambdaCentresIdx,:,j).*exp(-1i.*angle(sum(PSPsin(lambdaCentresIdx,:,j).*conj(PSPsin(lambdaCentresIdx-1,:,j)))));
                end
                PSPsout(lambdaCentresIdx,:,:) = Jc*squeeze(PSPsin(lambdaCentresIdx,:,:));
            end
        end
        
        [minV lambda0Idx] = min(abs(lambda0-lambdaShort));
        % fprintf('Centre Wavelength = %3.3f\n',lambda(lambda0Idx));
        [minV lambdaPlusIdx] = min(abs((lambda0+dLambda)-lambdaShort));
        [minV lambdaMinusIdx] = min(abs((lambda0-dLambda)-lambdaShort));
        dL = abs(lambda0Idx-lambdaPlusIdx);
        
        [minV lambda0Idx] = min(abs(lambda0-lambda));

        %Correlation for each basis
        ResponseSpectral_LG = zeros(1,totalModeCount,'single');
        ResponseSpectral_PS = zeros(lambdaShortCount,totalModeCount,'single');
        %Fourier transform of the correlation for each basis
        ResponseTemporal_PS = zeros(lambdaShortCount,totalModeCount,'single');
        ResponseTemporal_LG = zeros(lambdaShortCount,totalModeCount,'single');
        
        %Centre Wavelength transfer matrix
        Jc  = squeeze(MTMs(lambda0Idx,:,:));
        lambdaIdx=0;
        for lambdaIDX=1:lambdaStep:lambdaCount
            J = squeeze(MTMs(lambdaIDX,:,:));
            lambdaIdx=lambdaIdx+1;
            for modeIn=1:totalModeCount
                Pout = J*PSPin(:,modeIn);
                Pout = Pout./sqrt(sum(sum(abs(Pout).^2)));
                Jin = J(:,modeIn);
                Jin = Jin./sqrt(sum(sum(abs(Jin).^2)));
                Jcin = Jc(:,modeIn);
                Jcin = Jcin./sqrt(sum(sum(abs(Jcin).^2)));
                ResponseSpectral_LG(lambdaIdx,modeIn)     = sum((Jin.*conj(Jcin)));
                PSPoutNorm = squeeze(PSPout(:,modeIn));
                PSPoutNorm = PSPoutNorm./sqrt(sum(sum(abs(PSPoutNorm).^2)));
                ResponseSpectral_PS(lambdaIdx,modeIn)	= sum((Pout.*conj(PSPoutNorm)));
            end
        end

        %Apply a 250GHz Gaussian filter to the Fourier-transform of the
        %correlation function. Will be used to calculate the group-delay at
        %this wavelength.
        fc = c0./(lambda0.*1e-9);
        f = c0./(lambdaShort.*1e-9);defdef
        w = 250e9;
        gaussFilter = exp(-(f-fc).^2./(2.*w.^2));
        for modeIn=1:totalModeCount
            fD = squeeze(ResponseSpectral_PS(:,modeIn));
            tD = fftshift(fft(fftshift(fD.*gaussFilter)));
            ResponseTemporal_PS(:,modeIn) = tD;
            
            fD = squeeze(ResponseSpectral_LG(:,modeIn));
            tD = fftshift(fft(fftshift(fD.*gaussFilter)));
            ResponseTemporal_LG(:,modeIn) = tD;
        end
        
        ResponseTemporal_PSabs2 = abs(ResponseTemporal_PS).^2;
        ResponseTemporal_PSabs2 = ResponseTemporal_PSabs2./max(max(ResponseTemporal_PSabs2));
        ResponseTemporal_LGabs2 = abs(ResponseTemporal_LG).^2;
        ResponseTemporal_LGabs2 = ResponseTemporal_LGabs2./max(max(ResponseTemporal_LGabs2));
        
        GD_CentreOfMass = zeros(1,totalModeCount);
        for modeIn=1:totalModeCount
            GD_CentreOfMass(modeIn) = sum((ResponseTemporal_PSabs2(:,modeIn).*timeShort.'))./sum(ResponseTemporal_PSabs2(:,modeIn));
        end
        
        %Calculate the delta angular frequency required to calculate group
        %delays
        GD_dw = zeros(1,lambdaShortCount);
        for lambdaIdx=1:lambdaShortCount
            for modeIn=1:totalModeCount
                Pidx0 = lambdaIdx+dL;       
                Pidx1 = lambdaIdx-dL;
                if (Pidx0<=lambdaShortCount && Pidx1>=1)                    
                    f0 = c0./(lambdaShort(Pidx0).*1e-9);
                    f1 = c0./(lambdaShort(Pidx1).*1e-9);
                    dw = (2.*pi.*(f1-f0));
                    GD_dw(lambdaIdx) = dw;
                end
            end
        end

        [minV idx] = (min(abs(lambdaShort-lambda0)));

        GD = angle(GD);
        GroupDelay_Eigenvalues(lambdaCentresIdx,:) = GD;
        GroupDelay_Eigenvalues_dw(lambdaCentresIdx) = GD_dw(idx);
        GroupDelay_Eigenvalues_ps = GroupDelay_Eigenvalues;
        s = size(GroupDelay_Eigenvalues_ps);
        GroupDelay_CentreOfMass(lambdaCentresIdx,:) = GD_CentreOfMass;
        gdIdx = 1;
        if (lambdaCentresIdx==gdIdx)
            GroupDelay_CentreOfMassx = GroupDelay_CentreOfMass(gdIdx,:);
            %   GroupDelay_Eigenvalues_dw
            %  for i=1:s(2)
            GroupDelay_CentreOfMassx(gdIdx,:) = (GroupDelay_CentreOfMass(gdIdx,:).*GroupDelay_Eigenvalues_dw(gdIdx));
            % end
            dGD23 = round((GroupDelay_CentreOfMassx-GroupDelay_Eigenvalues(gdIdx,:))./(2.*pi));
            for i=1:s(2)
                GroupDelay_Eigenvalues(:,i) = GroupDelay_Eigenvalues(:,i)+2.*pi.*dGD23(i);
            end
        end
        for i=1:s(2)
            GroupDelay_Eigenvalues_ps(:,i) = (unwrap(GroupDelay_Eigenvalues(:,i))./GroupDelay_Eigenvalues_dw).*1e12;
        end
        if (lambdaCentresIdx==1)
            [Y sortIdx] = sort(GD_CentreOfMass);
            sortIdx=sortIdx(end:-1:1)
            GroupDelay_CentreOfMass = GroupDelay_CentreOfMass(:,sortIdx);
            GroupDelay_Eigenvalues_ps = GroupDelay_Eigenvalues_ps(:,sortIdx);
%            GDs = GDs(:,sortIdx);
            GD = GD(sortIdx);
            PSPout = PSPout(:,sortIdx);
            PSPin = PSPin(:,sortIdx);
            PSPsout(lambdaCentresIdx,:,:) = PSPsout(lambdaCentresIdx,:,sortIdx);
            PSPsin(lambdaCentresIdx,:,:) = PSPsin(lambdaCentresIdx,:,sortIdx);
            ResponseTemporal_PSabs2 = ResponseTemporal_PSabs2(:,sortIdx);
            GroupDelay_Eigenvalues = GroupDelay_Eigenvalues(:,sortIdx);
            ResponseSpectral_PS = ResponseSpectral_PS(:,sortIdx);
        end
        if ((max(max(GroupDelay_CentreOfMass)).*TUNITS)>TLIM(2))
            TLIM(2) = (max(max(GroupDelay_CentreOfMass))).*TUNITS+temporalWindowSize;
        end
        if ((min(min(GroupDelay_CentreOfMass)).*TUNITS)<TLIM(1))
            TLIM(1) = min(min(GroupDelay_CentreOfMass)).*TUNITS;
            TLIM(1) = TLIM(1)-temporalWindowSize;
         end   
        %TLIM = [-90 -70];
        figure(102021);
        subplot(1,2,1);
        plot(timeShort.*TUNITS,ResponseTemporal_PSabs2,'LineWidth',lineWidth);
        xlim(TLIM);

        [Y I] = sort(lambdaCentres);
        axis square;
        grid on;
        xlabel('time (ps)');
        ylabel('Response');
        set(gca,'LineWidth',lineWidth');
                        subplot(1,2,2);
        plot(timeShort.*TUNITS,10.*log10(ResponseTemporal_PSabs2),'LineWidth',lineWidth);
        xlim(TLIM);
        ylim([-40 0]);
        if (outputToVideo)
            fName = sprintf('%sResponseTime_%4.4i.png',videoOutputFolder,lambdaCentresIdx);
            print(fName,sprintf('-r%i',videoRes),'-dpng');
            
            if (lambdaCentresIdx==1)
                videoObjResponseTime = VideoWriter(fName,'Uncompressed AVI');
                videoObjResponseTime.FrameRate = fps;
                open(videoObjResponseTime);
            end
            fImg = imread(fName);
            f = im2frame(fImg);
            writeVideo(videoObjResponseTime,f);
        end
        figure(2);
       % subplot(1,4,1);
%        plot(lambdaShort,GDs.*1e12,'-x');
      %  axis square;
      %  title('Group Delay');
       % xlim(XLIM);
        
        ylabel('ps');
        xlabel('\lambda (nm)');
        axis square;
        grid on;
        subplot(1,4,2);
        % plot(1:totalModeCount,real(GD).^2,'-o',1:totalModeCount,imag(GD).^2,'-x',1:totalModeCount,abs(GD).^2,'-*');
        
        plot(1:totalModeCount,GD,'-o');
        axis square
        title('Group Delay');
        % ylim([0 1.1]);
        
        subplot(1,4,3);
        
        
        [Y I] = sort(lambdaCentres);
        plot(lambdaCentres(I),GroupDelay_Eigenvalues_ps(I,:),'-x');
        axis square;
        ylabel('ps');
        xlabel('\lambda (nm)');
        
        
        subplot(1,4,4);
        plot(lambdaCentres(I),GroupDelay_CentreOfMass(I,:).*TUNITS,'-x');
        axis square;
        
        GroupDelay_Eigenvaluesnorm = GroupDelay_Eigenvalues_ps(I,:);
        GroupDelay_CentreOfMassnorm = GroupDelay_CentreOfMass(I,:).*TUNITS;
        for modeIn=1:totalModeCount
            GroupDelay_Eigenvaluesnorm(:,modeIn) = GroupDelay_Eigenvaluesnorm(:,modeIn)-GroupDelay_Eigenvalues_ps(I,1);
            GroupDelay_CentreOfMassnorm(:,modeIn) = GroupDelay_CentreOfMassnorm(:,modeIn)-GroupDelay_CentreOfMass(I,1).*TUNITS;
        end
        figure(23423423);
        hold off;
        plot(lambdaCentres(I),GroupDelay_Eigenvaluesnorm,'x','LineWidth',lineWidth);
        hold on;
        plot(lambdaCentres(I),GroupDelay_CentreOfMassnorm,'-','LineWidth',lineWidth);
        xlim(XLIM);
        ylim(TLIM);
        if (TLIM(1)==-4 && TLIM(2)==4)
            set(gca,'YTick',-3:3);
        end
        set(gca,'XTick',lambdaTicks);
        set(gca,'XTickLabel',lambdaLabels);
        axis square;
        xlabel('\lambda (nm)');
        ylabel('Differential Delay (ps)');
        set(gca,'LineWidth',lineWidth');
        grid on;
        
        if (outputToVideo)
            fName = sprintf('%sGroupDelay_%4.4i.png',videoOutputFolder,lambdaCentresIdx);
            print(fName,sprintf('-r%i',videoRes),'-dpng');
            
            if (lambdaCentresIdx==1)
                videoObjGroupDelay = VideoWriter(fName,'Uncompressed AVI');
                videoObjGroupDelay.FrameRate = fps;
                open(videoObjGroupDelay);
            end
            
            fImg = imread(fName);
            f = im2frame(fImg);
            writeVideo(videoObjGroupDelay,f);
        end
        for modeIn=1:totalModeCount
            ResponseSpectral_LG(:,modeIn) = ResponseSpectral_LG(:,modeIn)./max(max(abs(ResponseSpectral_LG(:,modeIn))));
            ResponseSpectral_PS(:,modeIn) = ResponseSpectral_PS(:,modeIn)./max(max(abs(ResponseSpectral_PS(:,modeIn))));
        end
        
        
        h = figure(4);
        set(gcf,'DefaultAxesColorOrder',lineColor)
        ODB = 10.*log10(abs(ResponseSpectral_LG(:,:)).^2);
        ODB(:,ignoreModes) = nan;
        subplot(1,2,1);
        hold off;
        plot(lambdaShort,ODB,'LineWidth',lineWidth);
        hold on;
        plot([lambda0 lambda0],[-100 100],'-black','LineWidth',lineWidth);
        ylim(YLIM);
        xlim(XLIM);
        title('Laguerre-Gaussian');
        axis square;
        set(gca,'YTick',yTicks);
        set(gca,'XTick',lambdaTicks);
        set(gca,'XTickLabel',lambdaLabels);
        grid on;
        xlabel('\lambda (nm)');
        set(gca,'LineWidth',lineWidth);
        subplot(1,2,2);
        ResponseSpectral_PSDB = 10.*log10(abs(ResponseSpectral_PS).^2);
        ResponseSpectral_PSDB(:,ignoreModes) = nan;
        hold off;
        plot(lambdaShort,ResponseSpectral_PSDB,'LineWidth',lineWidth);
        hold on;
        plot([lambda0 lambda0],[-100 100],'-black','LineWidth',lineWidth);
        ylim(YLIM);
        worstDB = min(min(ResponseSpectral_PSDB));
        axis square;
        set(gca,'YTick',yTicks);
        set(gca,'XTick',lambdaTicks);
        set(gca,'XTickLabel',lambdaLabels);
        grid on;
        xlabel('\lambda (nm)');
        set(gca,'LineWidth',lineWidth);
        title(sprintf('Principal States'));
        xlim(XLIM);
        
        if (outputToVideo)
            fName = sprintf('%sResponseSpectral_%4.4i.png',videoOutputFolder,lambdaCentresIdx);
            print(fName,sprintf('-r%i',videoRes),'-dpng');
            
            if (lambdaCentresIdx==1)
                videoObjResponseSpectral = VideoWriter(fName,'Uncompressed AVI');
                videoObjResponseSpectral.FrameRate = fps;
                open(videoObjResponseSpectral);
            end
            
            fImg = imread(fName);
            f = im2frame(fImg);
            writeVideo(videoObjResponseSpectral,f);
        end
        
        
        spatialCoherence(lambdaCentresIdx,dLambdaIdx) = sum(sum(abs(ResponseSpectral_PS).^2))./sum(sum((ones(size(ResponseSpectral_PS)))));
        spatialCoherenceStates(lambdaCentresIdx,dLambdaIdx,:) = sum(abs(ResponseSpectral_PS).^2)./lambdaShortCount;
        spatialCoherenceWorst(lambdaCentresIdx,dLambdaIdx) = min(sum(abs(ResponseSpectral_PS).^2))./lambdaShortCount;
        spatialCoherenceLGStates(lambdaCentresIdx,dLambdaIdx,:) = sum(abs(ResponseSpectral_LG).^2)./lambdaShortCount;
        figure(24323432);
        stem(squeeze(spatialCoherenceStates(lambdaCentresIdx,dLambdaIdx,:)));
        ylim([0 1]);
        xlim([0 totalModeCount+1]);
        xlabel('Principal state');
        ylabel('Coherence');
        
        figure(24323439);
        hold off;
        plot(squeeze(spatialCoherenceStates(lambdaCentresIdx,dLambdaIdx,:)),'MarkerFaceColor',[0 0.5 0],'Marker','o','LineStyle','none','Color',[0 0.5 0]);
        hold on;
        plot(squeeze(spatialCoherenceLGStates(lambdaCentresIdx,dLambdaIdx,:)),'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none','Color',[1 0 0]);
        hold off;
        ylim([0 1]);
        xlim([0 totalModeCount+1]);
        xlabel('Principal state');
        ylabel('Coherence');
        axis square;
        grid on;
        
        if (lambdaCentresCount>1)
            figure(24323433);
            imagesc(1:totalModeCount,lambdaCentres,squeeze(spatialCoherenceStates));
        xlabel('Principal state');
        ylabel('Coherence');
            set(gca,'CLim',[0 1]);
            colormap(CMP_JET);
        end


        if (length(dLAMBDA)>1)
            figure(101);
            plot(dLAMBDA,spatialCoherence,'-x',dLAMBDA,spatialCoherenceWorst,'-o');
            xlabel('\Delta\lambda');
            ylabel('Coherence');
            [maxV maxIdx] = max(spatialCoherence);
            title(sprintf('%3.3f\n',dLAMBDA(maxIdx)));
        end

        %Plots images of the principal states
        if (visualiseStates || outputToVideo)
            for i=1:totalModeCount
                maxCount = 12;
                if (totalModeCount>maxCount)
                    
                else
                    maxCount = totalModeCount;
                end
                iOffset = floor((i-1)/maxCount);
                if (visualiseStates)
                    %Plot Principal mode response
                    figure(301+iOffset);
                    subplot(3,maxCount,mod((i-1),maxCount)+1);
                    
                    plot(lambdaShort,squeeze(ResponseSpectral_PSDB(:,i)));
                    ylabel('PM');
                    axis square;
                    ylim(YLIM);
                    xlim(XLIM);
                    grid on;
                end
                c = squeeze(PSPout(:,i));
                H = 0;
                V = 0;
                HV = 0;
                k=1;
                for j=1:2:totalModeCount
                    H = H+c(j).*squeeze(MODES(k,:,:));
                    V = V+c(j+1).*squeeze(MODES(k,:,:));
                    k=k+1;
                end
                s = size(H);
                
                HV = zeros(2.*s(1),s(2));
                HV(1:s(1),:) = H;
                HV((s(1)+1:end),:) = V;
                
                
                if (visualiseStates)

                    %Plot what principal mode looks like at output
                    subplot(3,maxCount,mod(i-1,maxCount)+1+2.*maxCount);
                    imagesc(abs(HV).^2);
                    title(sprintf('b_%i',i));
                    axis off;
                    axis equal
                end
                HV = abs(HV).^2;
                
                if (outputToVideo)
                    HV = HV./max(max(HV));

                    HV = uint8(255.*HV);
                    f = im2frame(HV,CMP_JET);

                    writeVideo(videoObjsOut{i},f);
                end
                
                c = squeeze(PSPin(:,i));
                H = 0;
                V = 0;
                HVin = 0;
                k=1;
                for j=1:2:totalModeCount
                    H = H+c(j).*squeeze(MODES(k,:,:));
                    V = V+c(j+1).*squeeze(MODES(k,:,:));
                    k=k+1;
                end
                s = size(H);
                
                HVin = zeros(2.*s(1),s(2));
                HVin(1:s(1),:) = H;
                HVin((s(1)+1:end),:) = V;
                
                if (visualiseStates)
                    %Principal mode at input
                    subplot(3,maxCount,mod(i-1,maxCount)+1+maxCount);
                    imagesc(abs(HVin).^2);
                    axis off;
                    axis equal
                    title(sprintf('a_%i',i));
                end
                % HVinFrame = abs(HVin).^2;
                HVin = abs(HVin).^2;
                
                if (outputToVideo)
                    HVin = HVin./max(max(HVin));
                    HVin = uint8(255.*HVin);
                    f = im2frame(HVin,CMP_JET);

                    writeVideo(videoObjsIn{i},f);
                end
                
                if (visualiseStates)
                    figure(401+iOffset);
                    %LG response
                    subplot(3,maxCount,i);
                    plot(lambdaShort,squeeze(ODB(:,i)));
                    axis square;
                    ylim(YLIM);
                    xlim(XLIM);
                    grid on;
                    ylabel('LG');
                end
                H = 0;
                V = 0;
                HV = 0;
                idx = floor((i-1)/2)+1;
                
                if (mod(i-1,2))
                    H = H+squeeze(MODES(idx,:,:));
                    s = size(H);
                else
                    V = V+squeeze(MODES(idx,:,:));
                    s = size(V);
                end
                c = MTMs(lambdaCentresIdx0,:,i);
                HV = zeros(2.*s(1),s(2));
                HV(1:s(1),:) = H;
                HV((s(1)+1:end),:) = V;
                
                if (visualiseStates)
                    %LG mode input;
                    subplot(3,maxCount,mod(i-1,maxCount)+1+maxCount);
                    imagesc(abs(HV).^2);
                    axis off;
                    axis equal
                end
                H = 0;
                V = 0;
                HVout = 0;
                k=1;
                for j=1:2:totalModeCount
                    H = H+c(j).*squeeze(MODES(k,:,:));
                    V = V+c(j+1).*squeeze(MODES(k,:,:));
                    k=k+1;
                end
                s = size(H);
                HVout = zeros(2.*s(1),s(2));
                HVout(1:s(1),:) = H;
                HVout((s(1)+1:end),:) = V;
                
                if (visualiseStates)
                    %LG mode output
                    subplot(3,maxCount,mod(i-1,maxCount)+1+2.*maxCount);
                    imagesc(abs(HVout).^2);
                    axis off;
                    axis equal
                end
            end
        end

        
        SecondOrderCorrelation_PS = zeros(totalModeCount,length(lambdaCentres));
        
        for lambda2IDX=1:lambdaCentresCount
            for i=1:totalModeCount
                SecondOrderCorrelation_PS(i,lambda2IDX) = (sum(PSPsout(lambda2IDX,:,i).*conj(PSPsout(lambdaCentresIdx0,:,i))));
            end
        end
        
        [Y I] = sort(lambdaCentres);
        figure(SupFigOffset+4);
        subplot(1,3,3);
        plot(lambdaCentres(I),10.*log10(abs(SecondOrderCorrelation_PS(:,I).').^2),'-x');
        xlim([lambda(1) lambda(end)]);
        ylim([-20 1]);
        
    end
    
end


[Y I] = sort(lambdaCentres);

for i=1:length(lambdaCentres)
        figure(SupFigOffset+4);
        subplot(1,3,3);
    hold off;
    PDB = 10.*log10(abs(SecondOrderCorrelation_PS(:,I).').^2);
    PDB(:,ignoreModes) = nan;
    plot(lambdaCentres(I),PDB,'-x','LineWidth',lineWidth);
    hold on;
    plot(lambdaCentres(I(i)).*[1 1],[-100 100],'-black','LineWidth',lineWidth);
    xlim(XLIM);
    set(gca,'XTick',lambdaTicks);
    set(gca,'XTickLabel',lambdaLabels);
    ylim([-12 1]);
    axis square;
    grid on;
    xlabel('\lambda (nm)');
    ylabel('2nd-order Correlation (dB)');
    
    if (outputToVideo)
        fName = sprintf('%sSecondOrderCorrelation_%4.4i.png',videoOutputFolder,i);
        print(fName,sprintf('-r%i',videoRes),'-dpng');
        if (i==1)
            videoObjSecondOrderCorrelation = VideoWriter(fName,'Uncompressed AVI');
            videoObjSecondOrderCorrelation.FrameRate = fps;
            open(videoObjSecondOrderCorrelation);
        end
        
        fImg = imread(fName);
        f = im2frame(fImg);
        writeVideo(videoObjSecondOrderCorrelation,f);
    end
end

for lambdaCentresIdx=1:length(lambdaCentres)
    figure(SupFigOffset+4);
    subplot(1,3,1);
    hold off;
    plot(lambdaCentres(I),GroupDelay_Eigenvaluesnorm,'x','LineWidth',lineWidth);
    hold on;
    plot(lambdaCentres(I),GroupDelay_CentreOfMassnorm,'-','LineWidth',lineWidth);
    hold on;
    plot(lambdaCentres(I(lambdaCentresIdx)).*[1 1],[-100 100],'-black','LineWidth',lineWidth)
    xlim(XLIM);
    ylim(TLIM);
    %set(gca,'YTick',-3:3);
    set(gca,'XTick',lambdaTicks);
    set(gca,'XTickLabel',lambdaLabels);
    axis square;
    xlabel('\lambda (nm)');
    ylabel('Differential Delay (ps)');
    set(gca,'LineWidth',lineWidth');
    grid on;
    
    if (outputToVideo)
        fName = sprintf('%sGroupDelayFull_%4.4i.png',videoOutputFolder,lambdaCentresIdx);
        print(fName,sprintf('-r%i',videoRes),'-dpng');
        
        if (lambdaCentresIdx==1)
            videoObjGroupDelayFull = VideoWriter(fName,'Uncompressed AVI');
            videoObjGroupDelayFull.FrameRate = fps;
            open(videoObjGroupDelayFull);
        end
        
        fImg = imread(fName);
        f = im2frame(fImg);
        writeVideo(videoObjGroupDelayFull,f);
    end
end
if (outputToVideo)
    close (videoObjGroupDelayFull);
    
    for i=1:totalModeCount
        close(videoObjsIn{i});
    end
    for i=1:totalModeCount
        close(videoObjsOut{i});
    end
    
    close(videoObjSecondOrderCorrelation);
    close(videoObjResponseSpectral);
    close(videoObjResponseTime);
    close(videoObjGroupDelay);
end

figure(43224324);
plot(1e12.*GroupDelay_CentreOfMass,'MarkerFaceColor',[0 0.5 0],'Marker','o','LineStyle','none','Color',[0 0.5 0]);
set(gca,'XTick',MGtick);
ylim([-250 20])
xlim([0.5 totalModeCount+0.5]);
grid on;
%set(gcf, 'Position', [0 0 1800 600])

%Plot the average transfer matrix (Supplemental Figure 6)
figure(SupFigOffset+6);

subplot(1,2,1);
imagesc(squeeze(mean(abs(MTMs))));
axis square;
set(gca,'XTick',MGtick);
set(gca,'YTick',MGtick);
set(gca, 'XColor', 'r')
set(gca, 'YColor', 'r')
colormap(gray(256));
title('Average |U|');
grid on;
subplot(1,2,2);
UU=0;
for lambdaIdx=1:lambdaCount
    M0 = squeeze(MTMs(lambdaIdx,:,:));
    U2 = M0*M0';
    UU=UU+abs(U2);
end
imagesc(UU);
axis square;
set(gca,'XTick',MGtick);
set(gca,'YTick',MGtick);
set(gca, 'XColor', 'r')
set(gca, 'YColor', 'r')
title('Average|U*U^*|');
grid on;

if (doDOP==1)
    dopResponse;
end