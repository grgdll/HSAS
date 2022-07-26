import subprocess
import os
import glob
import numpy as np

def satconCall(satconpath, logfile, calfiles, outfolder, skip=False):
    filebase = logfile.split('\\')[-1].split('.')[0]
 
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
		
    outfiles = []
    print logfile
##    print calfiles
##    print outfolder
	
    for cf in calfiles:
##        print cf
	
        instname = cf.split('\\')[-1].split('.')[0]
##        print instname
        
        outfile = outfolder+'\\'+filebase+'_'+instname+'.dat'
        print outfile
        
        if not skip:
            result = subprocess.call([satconpath, logfile, outfile, cf, '-x','-tf','-n','-efc','-OW=y','-EW=y','-ds']) #

            if result==0:
                outfiles.append(outfile)
                
        else:
            outfiles.append(outfile)
            
    return outfiles




print "\n"
print "YOU NEED TO RUN THIS SCRIPT ON A WINDOW COMPUTER USING PYTHON FOR WINDOWS"
print "\n"

	
bdir = r'Z:\grg\AMT25\Data\Underway_Rrs'
bcal = r'Z:\Optics_group\Data\allcals4extraction\HyperSAS_config\PreAMT29' 
satconpath = r'C:\Program Files\Satlantic\SatCon\SatCon.exe'
calfiles =  [bcal+r'\HLD222.cal', bcal+r'\HSL222.cal', 
             bcal+r'\HLD223.cal', bcal+r'\HSL223.cal', 
             bcal+r'\HLD464.cal', bcal+r'\HSL464.cal', 
             bcal+r'\HED258.cal', bcal+r'\HSE258.cal']



dirin = np.sort(glob.glob(bdir + '\\2015*'))

for idoy in dirin[9:-1]:
    
    doy = os.path.basename(idoy)


    dirdoy = bdir + '\\'+doy
    

    dirout = dirdoy + '\\' + doy + '_Satcon_extracted_raw'

    fn = glob.glob(idoy + '\\' + doy + '_raw\\' + '\\2015*' )

    for logfile in fn:
        satconCall(satconpath, logfile, calfiles, dirout, skip=False)
        










