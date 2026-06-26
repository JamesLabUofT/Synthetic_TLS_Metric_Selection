import numpy as np
import math
from scipy.ndimage import convolve1d, maximum_filter

def computeLacunarity(vox):

  voxOG = np.array(vox)
  
  truncated = True
  
  k=max(voxOG.shape)
  L = np.zeros(k)
  for r in range(1,k+1):
    convoluted = voxOG.copy()
    for i in range(3):
      convoluted = convolve1d(convoluted, np.ones(r), axis = i,mode ='constant')
  
    #if the data was truncated to a specific size, we want to strip the convoluted results obtained by windows that went "over the edge" of the data. Think of a 1D array with 4 elements 1, 2, 3, 4. To process the convolution, we've padded arount the array such that it has a values of 0 on either sides of the 4 elements. 
    #If we then convolve with a 1D array of three elements 1, 1, 1, the result without stripping will be 3, 6, 9, 7. Notice that values 3 and 7 on the edges don't represent our data: they result from a 0 that was artificially placed there for the purposes of the calculation. 
  # if, on the other hand, the data was NOT truncated, then we can expect that there would be 0s on either sides of the array anyways, so we don't need to remove those convolution results that included the added 0s. I can't really think of an example where the data would not be truncated in SOME way. 
  # Note: in the case of a 1D array, we could have just used a 'valid' convolutation. However, since we are working in 3 dimensions, this would mean that the biggest window would be a cube with sides of size min(shape of data). More explanations below
  
    if truncated:
      x, y, z = convoluted.shape
      xlpad, xrpad = 0, 0
      ylpad, yrpad = 0, 0
      zlpad, zrpad = 0, 0
      
      # If we don't check if the box is greater than one of the dimensions, then we'll remove everything along that dimension
      # Which would remove everything. Since our data is not bounded by a cube (the length, height and depth are not 
      # necessarily equal), we still want to calculate the lacunarity with a window that has a size max(length, height, depth). Since the window IS a cube, that means we'll almost always end up with a window that is bigger than the data along at least one of the dimensions
      # We disregard the following for now: In the case of a truncated lidar scan, we assume that the sky is empty, and the x and y parts should be of the same size, therefore the added 0s don't add fake gaps, they are very real and represent the empty space over the trees.
      pad = math.floor(r/2)
      if r > x: #The window is bigger than the dimension, so we should only keep one value along that dimension
        # In particular, we want the value at the middle of the dimension since it will guarantee that the window will have captured all values
        # along that dimension
        xlpad, xrpad = x//2-1, x//2
      else: # The window is smaller than the dimension, so remove the convolution values that overlapped with
        #the padding
        if r % 2 == 0:
          xlpad, xlpad = r//2 - 1, -(r//2)
        else:
          xlpad, xlpad = pad, -(pad)
      if r > y:
        ylpad, yrpad = y//2-1, y//2
      else:
        if r % 2 == 0:
          ylpad, yrpad = r//2 - 1, -(r//2)
        else:
          ylpad, yrpad = pad, -pad
      if r > z:
        zlpad, zrpad = z//2-1, z//2
      else:  
        if r % 2 == 0:
          zlpad, zrpad = r//2 - 1, -(r//2)
        else:
          zlpad, zrpad = pad, -pad
      
      #array[:0] is empty whereas array[:None] = array, which is the behaviour that we want when the right pad = 0
      
      
      convoluted = convoluted[xlpad:None if xrpad==0 else xrpad, 
                              ylpad:None if yrpad==0 else yrpad,
                              zlpad:None if zrpad==0 else zrpad]
    
      
    M, nM = np.unique(convoluted, return_counts = True)
  
    QM = nM/convoluted.size
    
    Z1 = np.sum(M * QM)
    
    Z2 = np.sum(M**2 * QM)
  
    L[r-1] = Z2/(Z1**2)
  #Scott 2022:
  Rw = 1/k*np.sum(range(1, k+1) * L/L[0])
  h = 1 - 2*Rw/(1+k)
  
  #cooper 2022:
  LTotal = 1/L[0]*np.sum(L)
  return([h,LTotal])



