function [ ret ] = plotqcfigsshew( qcstruct, plottitle )
%plotqcfigs(

set(gcf, 'position', [10,10, 1000,800])
subplot(2,2,1)
fulltitle=strcat(plottitle, ' SFNR');
shewhart(qcstruct.SFNR_roi_nodrift, fulltitle, qcstruct.DateMatlab, 0, 1, 0);
subplot(2,2,2)
fulltitle=strcat(plottitle, ' SNR');
shewhart(qcstruct.SNR, fulltitle, qcstruct.DateMatlab, 0, 1, 0);
subplot(2,2,3)
fulltitle=strcat(plottitle, ' Drift');
shewhart(qcstruct.LinearDrift, fulltitle, qcstruct.DateMatlab, 0, 1, 0);
subplot(2,2,4)
fulltitle=strcat(plottitle, 'Max Slice Sig Change');
shewhart(qcstruct.SliceMaxSigChange, fulltitle, qcstruct.DateMatlab, 0, 1, 0);

ret=0;

end

