
#############################################################
# pippi: parse it, plot it
# ------------------------
# Program for creating plotting scripts for pippi.
#
# Author: Pat Scott (patscott@physics.mcgill.ca)
# Originally developed: March 2012
#############################################################

from pippi_utils import *
import subprocess

# Define script-specific pip file entries
doComparison = dataObject('plot_comparison',boolean)
postMeanOnPost = dataObject('plot_posterior_mean_on_posterior_pdf',boolean)
postMeanOnProf = dataObject('plot_posterior_mean_on_profile_like',boolean)
bestFitOnPost = dataObject('plot_best_fit_on_posterior_pdf',boolean)
bestFitOnProf = dataObject('plot_best_fit_on_profile_like',boolean)
doLegend1D = dataObject('legend_on_1D',int_list)
doLegend2D = dataObject('legend_on_2D',intuple_list)
legendLoc1D = dataObject('legend_locations_1D',string_dictionary)
legendLoc2D = dataObject('legend_locations_2D',int_pair_string_dictionary)
doColourbar = dataObject('plot_colourbar_2D',intuple_list)
doHistograms = dataObject('plot_as_histograms_1D',boolean)
legendLines = dataObject('extra_legend_lines',string_list)
plotSize = dataObject('plot_size',string)
blame = dataObject('blame',string)
colours = dataObject('colour_scheme',internal)
axisRanges = dataObject('axis_ranges',floatuple_dictionary)	
keys = keys+[doComparison,postMeanOnPost,postMeanOnProf,bestFitOnPost,
        bestFitOnProf,doColourbar,doLegend1D,doLegend2D,legendLoc1D,legendLoc2D,
        doHistograms,legendLines,blame,colours,axisRanges]
# Define pip file entries to be read from savedkeys file
labels = dataObject('quantity_labels',string_dictionary)
dataRanges = dataObject('axis_ranges',floatuple_dictionary)

# Constants
blameFractionalVerticalOffset = 1.2e-2
likeColourbarString = 'Profile likelihood ratio $\Lambda=\mathcal{L}/\mathcal{L}_\mathrm{max}$'
postColourbarString = 'Relative probability $P/P_\mathrm{max}$'
defaultLegendLocation = 'bl'

def script(filename):
  # input: 	filename = the name of the pip file

  print

  # Parse pip file
  getIniData(filename,keys)

  # Make sure that comparison is turned off if comaprison filename is missing
  if doComparison.value and secChain.value is None:
    print '  Warning: comparison curves requested but no comparison file specified.\n  Skipping comparison...\n'
    doComparison.value = False

  # Strip extensions off chain filenames 
  baseFilename = re.sub(r"\..?.?.?$", '', mainChain.value)
  if doComparison.value: secFilename = re.sub(r"\..?.?.?$", '', secChain.value)
  
  # Retrieve labels and data ranges saved in earlier parsing run
  getIniData([baseFilename+'_savedkeys.pip'],[labels,dataRanges])

  # set colour scheme if it is undefined
  if colours.value is None: colours.value = basic
    
  # Create 1D plotting scripts
  if oneDplots.value is not None:

    # Loop over requested plots
    for plot in oneDplots.value:

      print '    Writing scripts for 1D plots of quantity ',plot
         
      # Set up filenames
      currentBase = baseFilename+'_'+str(plot)
      currentBaseMinimal = re.sub(r'.*/', '', currentBase)
      if doComparison.value:
        currentSec = secFilename+'_'+str(plot)
        currentSecMinimal = re.sub(r'.*/', '', currentSec)
 
      # Get plot limits
      xtrema = dictFallback(axisRanges,dataRanges,plot)

      # Make profile likelihood plotting scripts
      if doProfile.value:

        # Get contours
        if contours.value is not None: contourLevels = getContours(baseFilename,plot,'like')
        # Get best-fit point
        if bestFitOnProf.value: bestFit = getCentralVal(baseFilename,plot,'like')
        # Get posterior mean
        if postMeanOnProf.value: postMean = getCentralVal(baseFilename,plot,'post')

        # Open plotting shell script file for writing
        outfile = smart_open(currentBase+'_like1D.bsh','w')

        outfile.close

      # Make posterior pdf plotting scripts      
      if doPosterior.value:

        # Get contours
        if contours.value is not None: 
          mainContourLevels = getContours(baseFilename,plot,'post')
          if doComparison.value: secContourLevels = getContours(secFilename,plot,'post')
        # Get best-fit point
        if bestFitOnPost.value: bestFit = getCentralVal(baseFilename,plot,'like')
        # Get posterior mean
        if postMeanOnPost.value: postMean = getCentralVal(baseFilename,plot,'post')
          
        # Open plotting shell script file for writing
        outfile = smart_open(currentBase+'_post1D.bsh','w')

        outfile.close

      # Make profile-posterior comparison plotting scripts
      if doProfile.value and doPosterior.value:

        if contours.value is not None: 
          likeContourLevels = getContours(baseFilename,plot,'like')
          postContourLevels = getContours(baseFilename,plot,'post')
 
        # Open plotting shell script file for writing
        outfile = smart_open(baseFilename+'_'+str(plot)+'_combo1D.bsh','w')

        # Always plot both best fit and posterior mean on comparison plot

        outfile.close


  # Create 2D plotting scripts
  if twoDplots.value is not None:

    # Loop over requested plots
    for plot in twoDplots.value:

      print '    Writing scripts for 2D plots of quantities ',plot
      
      # Set up filenames
      currentBase = baseFilename+'_'+'_'.join([str(x) for x in plot])
      currentBaseMinimal = re.sub(r'.*/', '', currentBase)
      if doComparison.value:
        currentSec = secFilename+'_'+'_'.join([str(x) for x in plot])
        currentSecMinimal = re.sub(r'.*/', '', currentSec)
 
      # Get plot limits
      xtrema = dictFallback(axisRanges,dataRanges,plot[0])
      ytrema = dictFallback(axisRanges,dataRanges,plot[1])
      xRange = xtrema[1] - xtrema[0]
      yRange = ytrema[1] - ytrema[0]

      # Determine plot size
      if plotSize.value is None or plotSize.value is '':
        if doColourbar.value is not None and plot in doColourbar.value:
          plotSizeInternal = '12.5cm x 4in'
        else:
          plotSizeInternal = '11cm x 4in'
      else:
          plotSizeInternal = plotSize.value

      # Make profile likelihood plotting scripts
      if doProfile.value:

        # Get contours
        if contours.value is not None: 
          contourLevels = getContours(baseFilename,plot,'like')
        
        # Open plotting shell script file for writing 
        outfile = smart_open(currentBase+'_like2D.bsh','w')
        outfile.write('#!/usr/bin/env bash\n')
        outfile.write('# This plot script created by pippi '+pippiVersion+' on '+datetime.datetime.now().strftime('%c')+'\n')
        outfile.write('ctioga2\\\n')
        outfile.write('  --name '+currentBaseMinimal+'_like2D')
        outfile.write('  --plot-scale \'1.3\'\\\n')
        outfile.write('  --page-size \''+plotSizeInternal+'\'\\\n')
        if doColourbar.value is not None and plot in doColourbar.value:
          outfile.write('  --frame-margins 0.15,0.19,0.05,0.15\\\n')
        else:
          outfile.write('  --frame-margins 0.15,0.1,0.05,0.15\\\n')
        outfile.write('  --xrange '+str(xtrema[0])+':'+str(xtrema[1])+'\\\n')
        outfile.write('  --yrange '+str(ytrema[0])+':'+str(ytrema[1])+'\\\n')
        outfile.write('  --ylabel \''+labels.value[plot[1]]+'\'\\\n')
        outfile.write('  --xlabel \''+labels.value[plot[0]]+'\'\\\n')
        outfile.write('  --label-style x /scale 1.0 /shift 0.15 --label-style y /scale 1.0 /shift 0.15\\\n') 
        outfile.write('  --xyz-map\\\n')
        if doColourbar.value is not None and plot in doColourbar.value:
          outfile.write('  --new-zaxis zvalues /location right /bar_size \'0.5cm\'\\\n')
          outfile.write('  --label-style zvalues /angle 270 /shift 0.4\\\n')        
        outfile.write('  --plot '+currentBaseMinimal+'_like2D.ct2@1:2:3 ')
        if doColourbar.value is not None and plot in doColourbar.value: outfile.write('/zaxis zvalues ')
        outfile.write('/color-map \''+colours.value.colourMap(contourLevels,'like')+'\'\\\n')
        if doComparison.value:
          # Do everything for comparison chain
          if contours.value is not None:
            # Plot contours
            outfile.write('  --plot '+currentSecMinimal+'_like2D.ct2@1:2:3 /fill-transparency 1\\\n')
            for contour in contourLevels:
              outfile.write('  --draw-contour '+contour+' /color '+colours.value.comparisonProfContourColour2D+
                            ' /style '+colours.value.comparisonContourStyle+' /width '+colours.value.lineWidth2D+'\\\n')
          if bestFitOnProf.value and colours.value.comparisonBestFitMarker is not None: 
            # Get best-fit point and plot it
            bestFit = getCentralVal(secFilename,plot,'like')
            outfile.write('  --draw-marker '+str(bestFit[0])+','+str(bestFit[1])+' '+
                          colours.value.comparisonBestFitMarker+' /color \''+colours.value.comparisonBestFitColour+
                          '\' /scale '+str(colours.value.comparisonBestFitMarkerScale)+' \\\n')
          if postMeanOnProf.value and colours.value.comparisonPostMeanMarker is not None:
            # Get posterior mean and plot it
            postMean = getCentralVal(secFilename,plot,'post')
            outfile.write('  --draw-marker '+str(postMean[0])+','+str(postMean[1])+' '+
                          colours.value.comparisonPostMeanMarker+' /color \''+colours.value.comparisonPostMeanColour+
                          '\' /scale '+str(colours.value.comparisonPostMeanMarkerScale)+' \\\n')
        outfile.write('  --plot '+currentBaseMinimal+'_like2D.ct2@1:2:3 /fill-transparency 1\\\n')
        if contours.value is not None: 
          # Plot contours
          for contour in contourLevels:
            outfile.write('  --draw-contour '+contour+' /color '+colours.value.mainProfContourColour2D+
                          ' /style '+colours.value.mainContourStyle+' /width '+colours.value.lineWidth2D+'\\\n')
        if doLegend2D.value is not None and plot in doLegend2D.value:
          # Write legend
          try: 
            legendLocation = legendLoc2D.value[plot[0]][plot[1]]
          except KeyError:
            legendLocation = defaultLegendLocation
          outfile.write('  --legend-inside \''+legendLocation+'\' /scale 1.0 /dy 1.0\\\n')
          if legendLines.value is not None: 
            for x in legendLines.value: outfile.write('  --legend-line \''+x+'\'\\\n')
          outfile.write('  --legend-line \'Prof.~likelihood\'\\\n')                     
        if bestFitOnProf.value: 
          # Get best-fit point and plot it
          bestFit = getCentralVal(baseFilename,plot,'like')
          outfile.write('  --draw-marker '+str(bestFit[0])+','+str(bestFit[1])+' '+
                        colours.value.mainBestFitMarker+' /color \''+colours.value.mainBestFitColour+
                        '\' /scale '+str(colours.value.mainBestFitMarkerScale)+' \\\n')
        if postMeanOnProf.value: 
          # Get posterior mean and plot it
          postMean = getCentralVal(baseFilename,plot,'post')
          outfile.write('  --draw-marker '+str(postMean[0])+','+str(postMean[1])+' '+
                        colours.value.mainPostMeanMarker+' /color \''+colours.value.mainPostMeanColour+
                        '\' /scale '+str(colours.value.mainPostMeanMarkerScale)+' \\\n')
        if blame.value is not None:
          blameYCoordinate = str(blameFractionalVerticalOffset * yRange + ytrema[1])
          outfile.write('  --draw-text '+str(xtrema[1])+','+blameYCoordinate+' \''+blame.value+'\' /scale 0.5 /justification right\\\n')
        if doColourbar.value is not None and plot in doColourbar.value:
          # Do labelling for colourbar
          outfile.write('  --y2 --plot '+currentBaseMinimal+'_like2D.ct2@1:2:3 /fill-transparency 1\\\n') 
          outfile.write('  --axis-style y /decoration ticks --yrange '+str(ytrema[0])+':'+str(ytrema[1])+'\\\n')
          outfile.write('  --ylabel \''+likeColourbarString+'\' /shift 3.5 /angle 180 /scale 0.8\\\n')
        outfile.close
        subprocess.call('chmod +x '+currentBase+'_like2D.bsh', shell=True)

      # Make posterior pdf plotting scripts      
      if doPosterior.value:

        # Get contours
        if contours.value is not None: 
          mainContourLevels = getContours(baseFilename,plot,'post')
          if doComparison.value: secContourLevels = getContours(secFilename,plot,'post')
       
        # Open plotting shell script file for writing
        outfile = smart_open(currentBase+'_post2D.bsh','w')
        outfile.write('#!/usr/bin/env bash\n')
        outfile.write('# This plot script created by pippi '+pippiVersion+' on '+datetime.datetime.now().strftime('%c')+'\n')
        outfile.write('ctioga2\\\n')
        outfile.write('  --name '+currentBaseMinimal+'_post2D')
        outfile.write('  --plot-scale \'1.3\'\\\n')
        outfile.write('  --page-size \''+plotSizeInternal+'\'\\\n')
        if doColourbar.value is not None and plot in doColourbar.value:
          outfile.write('  --frame-margins 0.15,0.19,0.05,0.15\\\n')
        else:
          outfile.write('  --frame-margins 0.15,0.1,0.05,0.15\\\n')
        outfile.write('  --xrange '+str(xtrema[0])+':'+str(xtrema[1])+'\\\n')
        outfile.write('  --yrange '+str(ytrema[0])+':'+str(ytrema[1])+'\\\n')
        outfile.write('  --ylabel \''+labels.value[plot[1]]+'\'\\\n')
        outfile.write('  --xlabel \''+labels.value[plot[0]]+'\'\\\n')
        outfile.write('  --label-style x /scale 1.0 /shift 0.15 --label-style y /scale 1.0 /shift 0.15\\\n') 
        outfile.write('  --xyz-map\\\n')
        if doColourbar.value is not None and plot in doColourbar.value:
          outfile.write('  --new-zaxis zvalues /location right /bar_size \'0.5cm\'\\\n')
          outfile.write('  --label-style zvalues /angle 270 /shift 0.4\\\n')        
        outfile.write('  --plot '+currentBaseMinimal+'_post2D.ct2@1:2:3 ')
        if doColourbar.value is not None and plot in doColourbar.value: outfile.write('/zaxis zvalues ')
        outfile.write('/color-map \''+colours.value.colourMap(mainContourLevels,'post')+'\'\\\n')
        if doComparison.value:
          # Do everything for comparison chain
          if contours.value is not None:
            # Plot contours
            outfile.write('  --plot '+currentSecMinimal+'_post2D.ct2@1:2:3 /fill-transparency 1\\\n')
            for contour in secContourLevels:
              outfile.write('  --draw-contour '+contour+' /color '+colours.value.comparisonPostContourColour2D+
                            ' /style '+colours.value.comparisonContourStyle+' /width '+colours.value.lineWidth2D+'\\\n')
          if bestFitOnPost.value and colours.value.comparisonBestFitMarker is not None: 
            # Get best-fit point and plot it
            bestFit = getCentralVal(secFilename,plot,'like')
            outfile.write('  --draw-marker '+str(bestFit[0])+','+str(bestFit[1])+' '+
                          colours.value.comparisonBestFitMarker+' /color \''+colours.value.comparisonBestFitColour+
                          '\' /scale '+str(colours.value.comparisonBestFitMarkerScale)+' \\\n')
          if postMeanOnPost.value and colours.value.comparisonPostMeanMarker is not None:
            # Get posterior mean and plot it
            postMean = getCentralVal(secFilename,plot,'post')
            outfile.write('  --draw-marker '+str(postMean[0])+','+str(postMean[1])+' '+
                          colours.value.comparisonPostMeanMarker+' /color \''+colours.value.comparisonPostMeanColour+
                          '\' /scale '+str(colours.value.comparisonPostMeanMarkerScale)+' \\\n')
        outfile.write('  --plot '+currentBaseMinimal+'_post2D.ct2@1:2:3 /fill-transparency 1\\\n')
        if contours.value is not None: 
          # Plot contours
          for contour in mainContourLevels: 
            outfile.write('  --draw-contour '+contour+' /color '+colours.value.mainPostContourColour2D+
                          ' /style '+colours.value.mainContourStyle+' /width '+colours.value.lineWidth2D+'\\\n')
        if doLegend2D.value is not None and plot in doLegend2D.value:
          # Write legend
          try: 
            legendLocation = legendLoc2D.value[plot[0]][plot[1]]
          except KeyError:
            legendLocation = defaultLegendLocation
          outfile.write('  --legend-inside \''+legendLocation+'\' /scale 1.0 /dy 1.0\\\n')
          if legendLines.value is not None: 
            for x in legendLines.value: outfile.write('  --legend-line \''+x+'\'\\\n')            
          outfile.write('  --legend-line \'Marg.~posterior\'\\\n')                     
        if bestFitOnPost.value: 
          # Get best-fit point and plot it
          bestFit = getCentralVal(baseFilename,plot,'like')
          outfile.write('  --draw-marker '+str(bestFit[0])+','+str(bestFit[1])+' '+
                        colours.value.mainBestFitMarker+' /color \''+colours.value.mainBestFitColour+
                        '\' /scale '+str(colours.value.mainBestFitMarkerScale)+' \\\n')
        if postMeanOnPost.value: 
          # Get posterior mean and plot it
          postMean = getCentralVal(baseFilename,plot,'post')
          outfile.write('  --draw-marker '+str(postMean[0])+','+str(postMean[1])+' '+
                        colours.value.mainPostMeanMarker+' /color \''+colours.value.mainPostMeanColour+
                        '\' /scale '+str(colours.value.mainPostMeanMarkerScale)+' \\\n')
        if blame.value is not None:
          blameYCoordinate = str(blameFractionalVerticalOffset * yRange + ytrema[1])
          outfile.write('  --draw-text '+str(xtrema[1])+','+blameYCoordinate+' \''+blame.value+'\' /scale 0.5 /justification right\\\n')
        if doColourbar.value is not None and plot in doColourbar.value:
          # Do labelling for colourbar
          outfile.write('  --y2 --plot '+currentBaseMinimal+'_post2D.ct2@1:2:3 /fill-transparency 1\\\n') 
          outfile.write('  --axis-style y /decoration ticks --yrange '+str(ytrema[0])+':'+str(ytrema[1])+'\\\n')
          outfile.write('  --ylabel \''+postColourbarString+'\' /shift 3.5 /angle 180 /scale 0.8\\\n')
        outfile.close
        subprocess.call('chmod +x '+currentBase+'_post2D.bsh', shell=True)
  
      # Make profile-posterior comparison plotting scripts
      if doProfile.value and doPosterior.value:

        if contours.value is not None: 
          likeContourLevels = getContours(baseFilename,plot,'like')
          postContourLevels = getContours(baseFilename,plot,'post')

        # Open plotting shell script file for writing
        outfile = smart_open(baseFilename+'_'+'_'.join([str(x) for x in plot])+'_combo2D.bsh','w')

        
        # Always plot best fit and posterior mean on comparison plot

        outfile.close

def getContours(baseFilename,plot,statistic): 
  # Construct dimensionality of plot and string indicating specific plot (if any)
  if type(plot) == list: 
    [dim, plot] = [str(len(plot)), '' if statistic == 'like' else '_'+'_'.join([str(x) for x in plot])]
  else:
    [dim, plot] = ['1', '' if statistic == 'like' else '_'+str(plot)]
  # Open contour file
  contourfile = safe_open(baseFilename+plot+'_'+statistic+dim+'D.contours')
  # Read contents
  fileContents = contourfile.readline()
  while fileContents[0] == '#': fileContents = contourfile.readline()
  #Shut it
  contourfile.close
  levels = fileContents.split()
  return levels

def getCentralVal(baseFilename,plot,statistic): 
  # Find central value (either best fit or posterior mean) for requested plot
  # Open .best file
  bestfile = safe_open(baseFilename+'.best')
  # Read contents
  fileContents = bestfile.readline()
  while fileContents[0] == '#': fileContents = bestfile.readline()
  fileContents = bestfile.readlines()
  # Shut it
  bestfile.close
  if statistic == 'like':
    # Extract best fit
    point = fileContents[1].split()
  elif statistic == 'post':
    # Extract posterior pdf
    point = fileContents[3].split()
  else:
    # Never get here
    sys.exit('Error: unrecognised statistic in pippi_script.getCentralVal.\nQuitting...')
  # Choose the coordinates corresponding to the axes of the current plot
  if type(plot) == list: 
    coordinates = [point[x] for x in plot]
  else:
    coordinates = point[plot]
  return coordinates

def dictFallback(risky,safe,key):
  # Try to extract entry corresponding to key from risky dataObject, otherwise use safe dataObject
  try:
    return risky.value[key]
  except (KeyError, TypeError):
    return safe.value[key]


