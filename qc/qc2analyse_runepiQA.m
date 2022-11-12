% qc1analyse.runepiQA.m
proc_output_path = '/home/sapje1/data_sapje1/QC/RoutineQC/proc_output';
proc_nii_path = '/home/sapje1/data_sapje1/QC/RoutineQC/proc_temp/nii_proc';
transfer_path = '/home/sapje1/data_sapje1/QC/RoutineQC/summary/fortransfer';

clear res niifiles niifiles2 ii

[status, out] = system(['ls ' proc_nii_path '/*.nii']);

if(~status)
    niifiles2=strsplit(out, '\n');
    niifiles = niifiles2(~cellfun('isempty', niifiles2))';
    for ii=1:numel(niifiles)
        disp(niifiles{ii})
        tmp1 = epiQA_Siemens(niifiles{ii});
        res(ii) = tmp1;
    end
    % Place epiQA outputfiles in the standard QC data structure
    system(['cp -f ' proc_nii_path '//*.pdf ' transfer_path]);
    system(['cp -f ' proc_nii_path '//*.txt ' transfer_path]);
    system(['mv -f ' proc_nii_path '//*.pdf ' proc_output_path]);
    system(['mv -f ' proc_nii_path '//*.txt ' proc_output_path]);
else
    disp('Nothing to process.')
end


