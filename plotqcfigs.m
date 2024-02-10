function [ ret ] = plotqcfigs( qcstruct, plottitle )
%plotqcfigs(

set(gcf, 'position', [10,10, 1600,800])
subplot(2,2,1)
fulltitle=strcat(plottitle, ' SFNR')
plot(qcstruct.DateMatlab, qcstruct.SFNR_roi_nodrift )
title(fulltitle)

subplot(2,2,2)
fulltitle=strcat(plottitle, ' SNR')
plot(qcstruct.DateMatlab, qcstruct.SNR )
title(fulltitle)

subplot(2,2,3)
fulltitle=strcat(plottitle, ' Drift')
plot(qcstruct.DateMatlab, qcstruct.LinearDrift)
title(fulltitle)

subplot(2,2,4)
fulltitle=strcat(plottitle, 'Max Slice Sig Change')
plot(qcstruct.DateMatlab, qcstruct.SliceMaxSigChange)
title(fulltitle)

ret=0;

end

