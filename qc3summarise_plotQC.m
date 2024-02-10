summarypath = '/home/sapje1/data_sapje1/QC/RoutineQC/summary/';
transferpath = [summarypath 'fortransfer'];
InvisibleFig = 0;  %InvisibleFig=1 makes figures invisible, for bash script

qcfile = strcat(summarypath, '3TEQC_epi.txt');
qc3tetab = readtable(qcfile);
qc3te = table2struct(qc3tetab, 'ToScalar', 1);
qc3te.DateMatlab = datetime(qc3te.Date, 'InputFormat', 'yy_MM_dd');

qcfile = strcat(summarypath, '3TWQC_epi.txt');
qc3twtab = readtable(qcfile);
qc3tw = table2struct(qc3twtab, 'ToScalar', 1);
qc3tw.DateMatlab = datetime(qc3tw.Date, 'InputFormat', 'yy_MM_dd');

qcfile = strcat(summarypath, 'ConnectomQC_epi.txt');
qc3tmtab = readtable(qcfile);
qc3tm = table2struct(qc3tmtab, 'ToScalar', 1);
qc3tm.DateMatlab = datetime(qc3tm.Date, 'InputFormat', 'yy_MM_dd');

qcfile = strcat(summarypath, '7TQC_epi.txt');
qc7ttab = readtable(qcfile);
qc7t = table2struct(qc7ttab, 'ToScalar', 1);
qc7t.DateMatlab = datetime(qc7t.Date, 'InputFormat', 'yy_MM_dd');


%%% SHEWHART CHARTS %%%
figure(1)
if(InvisibleFig) == 1
    set(gcf,'Visible','off');
end
plotqcfigsshew(qc3te, '3TE');
export_fig([summarypath 'Shewhart3TE.pdf']);

figure(2)
if(InvisibleFig) == 1
    set(gcf,'Visible','off');
end
plotqcfigsshew(qc3tw, '3TW');
export_fig([summarypath 'Shewhart3TW.pdf']);

figure(3)
if(InvisibleFig) == 1
    set(gcf,'Visible','off');
end
plotqcfigsshew(qc3tm, 'Connectom');
export_fig([summarypath 'ShewhartConnectom.pdf']);

figure(4)
if(InvisibleFig) == 1
    set(gcf,'Visible','off');
end
plotqcfigsshew(qc7t, '7T');
export_fig([summarypath 'Shewhart7T.pdf']);

%%% UNPROCESSED PLOTS %%%
figure(5)
if(InvisibleFig) == 1
   set(gcf,'Visible','off');
end
set(gcf, 'Position', [ 10    10   600   800 ]);
subplot(2,1,1);

cmd = ['cp ' summarypath '*.pdf ' transferpath ];
system(cmd);

