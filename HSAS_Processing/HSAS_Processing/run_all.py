#---python 3.6 
import subprocess
import os
import time

result = subprocess.run(['octave','rad_functions/read_hsas_dir.m'],stdout=subprocess.PIPE)
dir = result.stdout.decode("utf-8") 
dir = dir[:-1]
print('Data DIR: ',dir)
sub_dir = [name for name in os.listdir(dir)]
print('All sub_DIR: ', sub_dir)


#os.system('octave calibrate_hsas.m')
for fd in sub_dir:
    print('# Starting for ',fd,'!')
    os.system('octave proc_L0_L1.m '+ fd)
    os.system('octave proc_L2.m '+ fd)
    os.system('octave amt_plot_L2.m '+ fd)
    time.sleep(2)
    print('# Done for ',fd,'!')

# alternative shortcuts
# ls ../../../../Processed_TPcorrection/L1/ | xargs -n1 -P10 octave -qf amt_plot_L2.m - short cut for plot functiom

