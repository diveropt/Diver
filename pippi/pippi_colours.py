
#############################################################
# pippi: parse it, plot it
# ------------------------
# Colour scheme module for pippi.  Add your own at the end of
# this file.
#
# Author: Pat Scott (patscott@physics.mcgill.ca)
# Originally developed: March 2012
#############################################################

import re

permittedSchemes = {}

class colourScheme:
  # Class for pippi plotting colour schemes

  # Name of colour scheme
  name = ''

  # Default values for colours
  mainPostColour1D = 'Blue'
  mainProfColour1D = 'Red'
  mainPostContourColour2D = 'Black'
  mainProfContourColour2D = 'Black'
  comparisonPostColour1D = 'Grey'
  comparisonProfcolour1D = 'Grey'
  comparisonPostContourColour2D = 'Grey'
  comparisonProfContourColour2D = 'Grey'
  baseProfColourMap = '#fff--#fff(contour2)--#f45(contour1)--#612'
  basePostColourMap = '#fff--#fff(contour2)--#88f(contour1)--#229'

  # Default values for 2D contour plotting styles
  mainContourStyle = 'Solid'
  comparisonContourStyle = 'Solid'
  lineWidth2D = '0.9'

  # Default markers and their colours
  referenceMarker = 'TimesOpen'
  referenceMarkerScale = 1.5
  referenceMarkerFillColour = 'Goldenrod'
  referenceMarkerStrokeColour = 'Grey'

  mainBestFitMarker = 'Star'
  mainBestFitMarkerScale = 1.5
  mainBestFitColour = '#300'

  mainPostMeanMarker = 'Bullet'
  mainPostMeanMarkerScale = 0.9
  mainPostMeanColour = '#004'

  comparisonBestFitMarker = 'Star'
  comparisonBestFitMarkerScale = 1.5
  comparisonBestFitColour = 'Grey'

  comparisonPostMeanMarker = 'Bullet'
  comparisonPostMeanMarkerScale = 0.9
  comparisonPostMeanColour = 'Grey'

  def __init__(self,name):
    global permittedSchemes
    name = name.lower()
    self.name = name
    if permittedSchemes is None:
      permittedSchemes = {name:self}
    else:
      permittedSchemes[name] = self

  def colourMap(self,contours,kind):
    #Construct colourmap from base colour map and contour levels
    if kind == 'post':
      localColourMap = self.basePostColourMap
    elif kind == 'like':
      localColourMap = self.baseProfColourMap
    else:
      sys.exit('    Error: unrecognised type of colourmap requested.\n    Quitting...\n')
    for i, contour in enumerate(contours):
      localColourMap = re.sub(r'contour'+str(i+1), contour, localColourMap) 
    return localColourMap
    
# basic colour scheme
basic = colourScheme('basic')

# iceCube colour scheme
iceCube = colourScheme('iceCube')
iceCube.baseProfColourMap = '#fff--#fff(contour2)--#292(contour1)--#f55(contour1)--#000'
iceCube.basePostColourMap = '#fff--#fff(contour2)--#29d(contour1)--#f55(contour1)--#000'
iceCube.mainBestFitColour = 'Black'
iceCube.mainPostMeanColour = 'Black'

# iceCube3sig colour scheme
iceCube = colourScheme('iceCube3sig')
iceCube.baseProfColourMap = '#fff--#fff(contour3)--#292(contour2)--#fff(contour2)--#929(contour1)--#f55(contour1)--#000'
iceCube.basePostColourMap = '#fff--#fff(contour3)--#29d(contour2)--#fff(contour2)--#929(contour1)--#f55(contour1)--#000'
iceCube.mainBestFitColour = 'Black'
iceCube.mainPostMeanColour = 'Black'

# SBClassic colour scheme
SBClassic = colourScheme('SBClassic')
SBClassic.baseProfColourMap = '#fff--#fff(contour2)--#2f2(contour1)--#f33(0.5)--#000'
SBClassic.basePostColourMap = '#fff--#fff(contour2)--#95d(contour1)--#f33(0.5)--#000'
SBClassic.mainBestFitColour = 'Black'
SBClassic.mainPostMeanColour = 'Black'

# nightOfTheAllanachs colour scheme
nightOfTheAllanachs = colourScheme('nightOfTheAllanachs')

