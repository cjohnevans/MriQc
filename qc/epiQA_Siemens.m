function [Results, SFNR_nii, epimasknii] = epiQA_Siemens(fname)
% function [Results, SFNR_nii, epimasknii] = epiQA_Siemens(fname)
%  1) From NIFTI fname.nii.gz calculate SFNR, SNR, Drift, slice-slice variation
%  2) Store output in Results struct
%  3) Store summary PDF

%CJE 151213.  Change to load_nii for better viewing, but keep consistency
%with load_nifti
%epi=load_nifti(fname);
% This preproc now done in setup files - so assume it's *.nii
%epi = load_nii_gz(fname);
epi = load_nii(fname);
epi.vol = double(flipdim(epi.img, 1));

PSCepi=diff(epi.vol,1,4);
episize=size(PSCepi);

xdim=episize(1);
ydim=episize(2);
zdim=episize(3);
nt=episize(4);

% threshold using fsl
cmd = ['/cubric/software/fsl/bin/fslstats ' fname ' -r' ];
[imgmin, imgmax]=system(cmd);
imgmax=str2num(imgmax);
imgmax=imgmax(2);

thresh = imgmax/10;

[aa, meanepi]=system(['/cubric/software/fsl/bin/fslstats ' fname ' -l ' num2str(thresh) ' -M']);
meanepi=str2num(meanepi);

PSCepi=100*(PSCepi./meanepi); % percent signal change
%PSCpercentile75=PSCepi(PSCepi, 75, 4);

meanvol=mean(epi.vol,4);
meanvol(:,:,zdim:64)=0;
meansag=flipdim(permute(meanvol, [3 2 1]),1);
meanrow=reshape(meansag, [xdim ydim*64]);
meansag_2d=cat(1, meanrow(:,(1:8*ydim)), meanrow(:,(1+8*ydim):(16*ydim)),  ...
      meanrow(:,(1+16*ydim):(24*ydim)), meanrow(:,(1+24*ydim):(32*ydim)), ...
      meanrow(:,(1+32*ydim):(40*ydim)), meanrow(:,(1+40*ydim):(48*ydim)), ...
      meanrow(:,(1+48*ydim):(56*ydim)), meanrow(:,(1+56*ydim):(64*ydim)) );


% Signal to Fluctuation Noise Ratio, SFNR (Glover) 
% NEED TO DEAL WITH LOW SIGNAL REGIONS IN SFNR CALC
noise=std(epi.vol,0,4);
noise(:,:,zdim:64)=0; %pad to 64 slices

SFNR = meanvol ./ noise;
%SFNRrow=reshape(SFNR, [xdim ydim*64]);
%SFNR_2d=cat(1, SFNRrow(:,(1:8*ydim)), SFNRrow(:,(1+8*ydim):(16*ydim)),  ...
%      SFNRrow(:,(1+16*ydim):(24*ydim)), SFNRrow(:,(1+24*ydim):(32*ydim)), ...
%      SFNRrow(:,(1+32*ydim):(40*ydim)), SFNRrow(:,(1+40*ydim):(48*ydim)), ...
%      SFNRrow(:,(1+48*ydim):(56*ydim)), SFNRrow(:,(1+56*ydim):(64*ydim)) );

% Display four slices: midslice-7 midslice-2 midslice+3 midslice+8
midslice = round(episize(3)/2);
dispslice = [midslice-7 midslice midslice+1 midslice+8 ];
SFNR_axdisp = cat(2, cat(1, rot90(SFNR(:,:,dispslice(1))), rot90(SFNR(:,:,dispslice(2)))), ...
    cat(1, rot90(SFNR(:,:,dispslice(3))), rot90(SFNR(:,:,dispslice(4)))));
%Display ROIs for SFNR calc
SFNR_axdisp(22:44, (64+22)) = 0;
SFNR_axdisp(22:44, (64+44)) = 0;
SFNR_axdisp(22, (64+22):(64+44)) = 0;
SFNR_axdisp(44, (64+22):(64+44)) = 0;
SFNR_axdisp((64+22):(64+44), 22) = 0;
SFNR_axdisp((64+22):(64+44), 44) = 0;
SFNR_axdisp((64+22), 22:44) = 0;
SFNR_axdisp((64+44), 22:44) = 0;

SFNR_midsl1 = rot90(SFNR(:,:,midslice));
SFNR_midsl1(22:44,22) = 0;
SFNR_midsl1(22:44,44) = 0;
SFNR_midsl1(22,22:44) = 0;
SFNR_midsl1(44,22:44) = 0;

SFNR_midsl2 = rot90(SFNR(:,:,midslice+1));
SFNR_midsl2(22:44,22) = 0;
SFNR_midsl2(22:44,44) = 0;
SFNR_midsl2(22,22:44) = 0;
SFNR_midsl2(44,22:44) = 0;


%do the same for sagital reslicing
SFNRsag=flipdim(permute(SFNR, [3 2 1]),1);
SFNRrow=reshape(SFNRsag, [xdim ydim*64]);
SFNRsag_2d=cat(1, SFNRrow(:,(1:8*ydim)), SFNRrow(:,(1+8*ydim):(16*ydim)),  ...
      SFNRrow(:,(1+16*ydim):(24*ydim)), SFNRrow(:,(1+24*ydim):(32*ydim)), ...
      SFNRrow(:,(1+32*ydim):(40*ydim)), SFNRrow(:,(1+40*ydim):(48*ydim)), ...
      SFNRrow(:,(1+48*ydim):(56*ydim)), SFNRrow(:,(1+56*ydim):(64*ydim)) );

% PSC stdev plot - better at picking up the rf.  Needs testing in
% script
% what does the stdev of PSC actually mean?
PSCstd=std(PSCepi,1,4);
PSCstd(:,:,zdim:64)=0;
PSCstd=flipdim(permute(PSCstd, [3 2 1]),1); %reformat to sag (rf
                                            %shows up better )

% ... not really CVrow (but saves on typing)	    
PSCrow=reshape(PSCstd, [xdim ydim*64]);   
PSC_2d=cat(1, PSCrow(:,(1:8*ydim)), PSCrow(:,(1+8*ydim):(16*ydim)),  ...
   PSCrow(:,(1+16*ydim):(24*ydim)), PSCrow(:,(1+24*ydim):(32*ydim)), ...
   PSCrow(:,(1+32*ydim):(40*ydim)), PSCrow(:,(1+40*ydim):(48*ydim)), ...
   PSCrow(:,(1+48*ydim):(56*ydim)), PSCrow(:,(1+56*ydim):(64*ydim)) );

PSC_std_MIP = max(PSCstd, [],3);
PSC_std_MIP = PSC_std_MIP(35:end,:);

%average PSC across each slice for all timepoints
PSCslice=squeeze(mean(mean(abs(PSCepi),1),2));

% Noise image (diff between vol 1 and 2)
noise_img = (epi.vol(:,:,(midslice:midslice+1),end-1) - epi.vol(:,:,(midslice:midslice+1),end));
drift_img = rot90(epi.vol(:,:,midslice,1) - epi.vol(:,:,midslice,end));
noise_img2(1:64,1:64,1:2) = 0;

for kk=1:2
    for ii=1:64
        for jj=1:64
            if(noise_img(ii,jj,kk) == 0)
                noise_img2(ii,jj,kk) = NaN;
            else
                noise_img2(ii,jj,kk) = noise_img(ii,jj,kk);
            end
        end
    end
end


%Summary Stats
% calc SFNR_Summary_DF later... need drift fix.
SFNR_Summary = mean(mean(mean((SFNR(23:43,23:43,(midslice:midslice+1)))))); %Generate from 21x21x2 voxel at centre
Signal_Summary = mean(mean(mean(meanvol(:,:,(midslice:midslice+1)))));
SNR_Summary = Signal_Summary / (std(reshape(noise_img, [numel(noise_img) 1])));
SlicePSC_Max = max(max(PSCslice));
DriftArtefact=100*std(reshape(drift_img, [numel(drift_img) 1]))/Signal_Summary;
%scale noise and drift images
noise_img = 100*noise_img ./ Signal_Summary;
noise_img_plot = rot90(noise_img(:,:,1));
drift_img = 100 * drift_img ./Signal_Summary;


% Signal Timecourse
Signal_ROI = epi.vol(23:43,23:43,(midslice:midslice+1), :);
Signal_ROI = squeeze(mean(mean(mean(Signal_ROI, 3), 1),2));
Signal_ROI = (100*Signal_ROI ./ Signal_ROI(1))-100; % as percentage
% Linear fit of drift
%volno=[1:300];
volno=[1:length(Signal_ROI)];
[fitvals] = polyfit(volno, Signal_ROI',1);
x1= fitvals(1); x0 = fitvals(2);
DriftFit = x1*volno + x0;
LinearDrift=x1*300; % in % over whole timeseries

%Need to recalculate this... without drift
epiROI = epi.vol(23:43,23:43,(midslice:midslice+1),:);
epiROImean = squeeze(mean(mean(mean(epiROI)))); %Generate from 21x21x2 voxel at centre%meanvol=mean(epi.vol,4);
% 2nd order fit
[order2poly] = polyfit(volno/max(volno), epiROImean',2);
order2fit = polyval(order2poly, volno/max(volno));
order2fix = polyval([order2poly(1) order2poly(2) 0 ] , volno/max(volno)); % don't correct zeroth order term for the sfnr calc.
%Need to recalculate this... without drift
epiROI = epi.vol(23:43,23:43,(midslice:midslice+1),:);
epiROImean = squeeze(mean(mean(mean(epiROI)))); %Generate from 21x21x2 voxel at centre%meanvol=mean(epi.vol,4);
% 2nd order fit
[order2poly] = polyfit(volno/max(volno), epiROImean',2);
order2fit = polyval(order2poly, volno/max(volno));
DriftMINMAX = 100*(max(order2fit) - min(order2fit))/order2fit(1);
%plot(volno, epiROImean, volno, order2fit);
driftfix = reshape(-order2fix, 1,1,1,max(volno));
driftfixroi = repmat(driftfix, [21 21 2]);
driftfixvol = repmat(driftfix, [64 64 zdim]);
epiROIdf = epiROI + driftfixroi;
epiROIdf_plot = 100* (squeeze(mean(mean(mean(epiROIdf)))) ... 
    ./ squeeze(mean(mean(mean(epiROIdf(:,:,:,1)))))) - 100; % for plot, as %
SFNR_DF =  mean(epiROIdf, 4) ./ std( epiROIdf, [],4);
SFNR_Summary_DF = mean(mean(mean(SFNR_DF)));
epivolDfix = epi.vol + driftfixvol;
SFNRvolDfix = mean(epivolDfix,4) ./ std(epivolDfix, 0,4);

% Whole volume SFNR, mask at 10% of max signal
mask_lim = 0.1 * max(max(max(epi.vol(:,:,:,1))));
epimask=epi.vol(:,:,:,1);
epimask(epimask<mask_lim)=0;
epimask(epimask>0)=1;
epimasknii=make_nii(epimask);
SFNRvolDfix = SFNRvolDfix .* epimask;
SFNRallvals = reshape(SFNRvolDfix, [64*64*zdim 1]);
SFNRallvals(SFNRallvals==0) = [];
SFNRallvalmean = mean(SFNRallvals);
SFNRallvalsem = std(SFNRallvals) ./ (length(SFNRallvals).^0.5);
%figure(501);plot(SFNRallvals)

% SNR_vol calcs
noiseVol_img = (epi.vol(:,:,:,end-1) - epi.vol(:,:,:,end));
sigVol_img = 0.5 * (epi.vol(:,:,:,end-1) + epi.vol(:,:,:,end));
sigVol_img = sigVol_img .* epimask;
SNRallvals = reshape(sigVol_img, [64*64*zdim 1]);
SNRallvals(SNRallvals==0) = [];
SNRallvalmean = mean(SNRallvals);
SNR_vol =  SNRallvalmean  ./ std(reshape(noiseVol_img, [numel(noiseVol_img) 1]));


% SFNR from mid slices
h=figure;
set(h,'Visible','off');
%set(h, 'Position', [ 134    87   560   887 ]);
set(h, 'Position', [ 10    10   600   800 ]);
subplot(4,2,1)
imagesc(SFNR_midsl1)
axis square
txt = ['SFNR (Vol Scale)' ];
title(txt);
colorbar
axis off

 subplot(4,2,2)
 clim=[0 2*SFNR_Summary_DF];
imagesc(SFNR_midsl1, clim)
axis square
txt = ['SFNR (2*ROI Scale)' ];
title(txt);
colorbar
axis off

 subplot(4,2,4)
%imagesc(rot90(meanvol(:,:,midslice)));
% imagesc(rot90(epi.vol(:,:,midslice,1)));
% axis square
% txt = ['Typical Image (slice 1)' ];
% title(txt);
% colorbar
% axis off
plot(volno, Signal_ROI, 'b', volno, epiROIdf_plot, 'r');
txt = ['Signal Drift (%)' ];
title(txt);
xlabel('Volume');
%ylabel('Signal')
legend('Uncorrected','Corrected', 'Location', 'SouthWest');

subplot(4,2,7)
imagesc(PSCslice)
colorbar
colormap(hot)
xlabel('Volume')
title('Slice % Signal Change')
ylabel('Slice')

subplot(4,2,5)
imagesc(noise_img_plot)
axis square
colorbar
%colorbar('location','southoutside')
axis off
title('Noise Image')

subplot(4,2,6)
imagesc(drift_img)
colorbar
%colorbar('location','southoutside')
axis off
title('Drift Image')
axis square

subplot(4,2,3)
imagesc(PSC_std_MIP, [0 max(PSC_std_MIP(:))]);
colormap hot
axis square
title('Fluctuation Map (MIP)');
colorbar
axis off

subplot(4,2,8)
axis off

[fpath, fname_short, fileext] = fileparts(fname);
tmp = [ fname_short ];
tmp = regexprep(tmp, '_','-');
text(0,0.9, tmp, 'FontName', 'Courier');

tmp = [ 'SFNR+Drift(ROI): ' num2str(SFNR_Summary,'%.1f') ];
text(0,0.8, tmp, 'FontName', 'Courier');

tmp = [ 'SFNR(ROI,Vol)  : ' num2str(SFNR_Summary_DF, '%.1f') ', ' num2str(SFNRallvalmean, '%.1f') ];
text(0,0.7, tmp, 'FontName', 'Courier');

tmp = [ 'SNR(ROI, Vol)  : ' num2str(SNR_Summary,'%.1f') ', ' num2str(SNR_vol, '%.1f') ];
text(0,0.6, tmp, 'FontName', 'Courier');

tmp = [ 'Signal         : ' num2str(Signal_Summary, '%.1f') ];
text(0,0.5, tmp, 'FontName', 'Courier');

tmp = [ 'Max \Delta Slice    : ' num2str(SlicePSC_Max, '%.3f') ' %' ];
text(0,0.4, tmp, 'FontName', 'Courier');

tmp = [ 'Drift%         : ' num2str(DriftMINMAX, '%.3f') ];
text(0,0.3, tmp, 'FontName', 'Courier');


%pdfname=[ fname(1:(end-7)) '.pdf' ];
pdfname = [fpath '//' fname_short '.pdf'];
% Use export_fig from matlabcentral
export_fig(pdfname);
outd = ['EPI QA analysis saved to ' pdfname ];
disp(outd)

%work out filename, date and time.  Assumes filename of the format
% QA3TE-GloverGSQAP-20_03_16-14_15_48-STD-1_3_12_2_1107_5_2_43_66075
fname_short_split = strsplit(fname_short, '-');

Results.File = fname_short;
Results.Scanner = fname_short_split{1};
Results.Date = fname_short_split{3};
Results.Time = fname_short_split{4};
Results.SNR_roi = round(SNR_Summary, 2); %roi
Results.SNR_vol = round(SNR_vol,2);
Results.SFNR_roi = round(SFNR_Summary,2);
Results.SFNR_roi_nodrift = round(SFNR_Summary_DF,2);
Results.SFNR_vol_nodrift = round(SFNRallvalmean, 2);
Results.semSFNR_vol_nodrift = round(SFNRallvalsem,4);
Results.LinearDrift = round(LinearDrift, 4);
Results.SliceMaxSigChange = round(SlicePSC_Max,4);
Results.DriftArtefact = round(DriftArtefact, 4);
SFNR_nii = make_nii(SFNRvolDfix);

writetable(struct2table(Results), [fpath '//' fname_short '.txt']);

