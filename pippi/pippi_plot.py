
#############################################################
# pippi: parse it, plot it
# ------------------------
# Plotting program for pippi.
#
# Author: Pat Scott (patscott@physics.mcgill.ca)
# Originally developed: March 2012
#############################################################

from pippi_utils import *

#Define plot-specific pip file entries
outdir = dataObject('output_dir',safe_string)
prepend = dataObject('file_prefix',safe_string)
append = dataObject('file_suffix',safe_string)
keys = keys+[outdir,prepend,append]

def plot(filename):

  print

  # Parse pip file
  getIniData(filename,keys)

  # Split into directory and filename without extension
  baseFilename = [re.sub(r'/.*?$', '/', mainChain.value), 
                  re.sub(r'.*/|\..?.?.?$', '', mainChain.value)]     
  
  # Run 1D plotting scripts
  if oneDplots.value is not None:
    # Work through 1D plotting scripts
    for plot in oneDplots.value:
      print '    Running plotting scripts for 1D plots of quantity ',plot       
      # Set up filenames
      currentBase = baseFilename[1]+'_'+str(plot)
      # Make profile likelihood plots
      if doProfile.value: 
        subprocess.call('cd '+baseFilename[0]+'; ./'+currentBase+'_like1D.bsh', shell=True)
        subprocess.call('mv '+baseFilename[0]+currentBase+'_like1D.pdf '+
         outdir.value+'/'+prepend.value+currentBase+'_like1D'+append.value+'.pdf', shell=True)
      # Make posterior pdf plots      
      if doPosterior.value: 
        subprocess.call('cd '+baseFilename[0]+'; ./'+currentBase+'_post1D.bsh', shell=True)
        subprocess.call('mv '+baseFilename[0]+currentBase+'_post1D.pdf '+
         outdir.value+'/'+prepend.value+currentBase+'_post1D'+append.value+'.pdf', shell=True)
      # Make profile-posterior comparison plots
      #if doProfile.value and doPosterior.value:
        #subprocess.call('cd '+baseFilename[0]+'; ./'+currentBase+'_combo1D.bsh', shell=True)
        #subprocess.call('mv '+baseFilename[0]+currentBase+'_combo1D.pdf '+
        # outdir.value+'/'+prepend.value+currentBase+'_combo1D'+append.value+'.pdf', shell=True)

  # Create 2D plotting scripts
  if twoDplots.value is not None:
    # Loop over requested plots
    for plot in twoDplots.value:
      print '    Running plotting scripts for 2D plots of quantity ',plot    
      # Set up filenames
      currentBase = baseFilename[1]+'_'+'_'.join([str(x) for x in plot])      
      # Make profile likelihood plots
      if doProfile.value: 
        subprocess.call('cd '+baseFilename[0]+'; ./'+currentBase+'_like2D.bsh', shell=True)
        subprocess.call('mv '+baseFilename[0]+currentBase+'_like2D.pdf '+
         outdir.value+'/'+prepend.value+currentBase+'_like2D'+append.value+'.pdf', shell=True)
      # Make posterior pdf plots      
      if doPosterior.value: 
        subprocess.call('cd '+baseFilename[0]+'; ./'+currentBase+'_post2D.bsh', shell=True)
        subprocess.call('mv '+baseFilename[0]+currentBase+'_post2D.pdf '+
         outdir.value+'/'+prepend.value+currentBase+'_post2D'+append.value+'.pdf', shell=True)
      # Make profile-posterior comparison plots
      #if doProfile.value and doPosterior.value:
        #subprocess.call('cd '+baseFilename[0]+'; ./'+currentBase+'_combo2D.bsh', shell=True)
        #subprocess.call('mv '+baseFilename[0]+currentBase+'_combo2D.pdf '+
        # outdir.value+'/'+prepend.value+currentBase+'_combo2D'+append.value+'.pdf', shell=True)

