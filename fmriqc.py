#!/home/john/miniconda3/envs/mri/bin/python
import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os

def splash():
    print('\nfmriqc.pc')
    print('   Rudamentory fmri quality control check\n')
    print('   USAGE:        python path_to_fmriqc/fmriqc.py FMRI_NIFTI ISINVIVO\n')
    print('   ARGUMENTS:    FMRI_NIFTI - 4D data to check (nii or nii.gz)')
    print('                 ISINVIVO   - y/n  y if in vivo data, n if phantom\n')
    return

# argv[0] is command, argv[1] .. argv[n] are arguments
if len(sys.argv) < 3:
    splash()
    sys.exit()
if not os.path.exists(sys.argv[1]):
    print("Can't find file " + sys.argv[1])
    sys.exit()
else:
    filename = sys.argv[1]
    if sys.argv[2] == 'y' or sys.argv[2] == 'Y':
        in_vivo = True
    if sys.argv[2] == 'n' or sys.argv[2] == 'N':
        in_vivo = False
    else:
        splash()
        print("!!! ISINVIVO should be 'y' or 'n'")
    fmri = mriqc.FmriQc(filename, in_vivo, run_report=True)


