from pylab import *
from scipy.ndimage import measurements
import numpy

def indices_max_all(a):
    # Only for 1D array (prob easy to extend)
    # Can also make a version that finds the 10 largest clusters or so and finds their indices
    indices = []
    maxa = max(a)
    for i in range(len(a)):
        if a[i]==maxa:
            indices.append(i)
    return indices

def indices_Nlargest(a,N):
    # Only for 1D array (prob easy to extend)
    # Can also make a version that finds the 10 largest clusters or so and finds their indices
    indices = []
    values  = []
    b = copy(a)
    for j in range(N):
        maxa = max(b)
        if j>0 and maxa<=1: # We don't want to look at very small clusters. 1 is the only cluster that is possible for Lx=Ly=1.
            break
        for i in range(len(b)):
            if b[i]==maxa:
                indices.append(i)
                values.append(maxa)
                b[i]=0
    return indices, values
    
def indices_Nlargest_delshift(indexarray, valuearray,k):
    # Only for 1D array (prob easy to extend)
    # Can also make a version that finds the 10 largest clusters or so and finds their indices
    # k is the element we want to remove from both arrays
    #index = k
    removedelement = indexarray[k]
    outindexarray  = zeros(len(indexarray)-1)
    outvaluearray  = zeros(len(indexarray)-1)
    counter = 0
    for i in range(len(indexarray)):
        if i!=k:
            if indexarray[i]>removedelement: # If the index is bigger than a the index of a removed element, we need to adjust for the removal
                indexarray[i]-= 1
            outindexarray[counter] = indexarray[i]
            outvaluearray[counter] = valuearray[i]
            counter += 1
    return outindexarray, outvaluearray

def findP_av(Lx,Ly,p,N_real):
    P_av = 0
    p = float(p)
    for i in range(N_real):
        r = rand(Lx,Ly)
        z = r<p
        #print(z)
        lw, num = measurements.label(z)
        
        #print(lw)
        b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
        shuffle(b) # shuffle this array
        shuffledLw = b[lw] # replace all values with values from b
        
        # Calculate areas
        area = measurements.sum(z, lw, index=arange(lw.max() + 1))
        areaImg = area[lw]
        
        '''
        figure()
        im3 = imshow(areaImg, origin='lower', interpolation='nearest')
        title("Clusters by area")
        '''
        
        #print(areaImg)
        #print(area)
        #print("area printed")
        #print("Type, area:", type(area))
        
        # Bounding box
        # The cluster with the largest area is not always the one with the largest bounding box
        maxarea = areaImg.max()
        #indices_maxarea = indices_max_all(area)      # Finds the index of the cluster(s) that have the largest area
        indices_maxarea, values_maxarea = indices_Nlargest(area,10)   # Finds the index of the clusters that have areas corresponding to the N largest occurances

        # Just assuming that the largest cluster is the spanning cluster: (Seems fair from comparing with the other method)
        if(len(indices_maxarea)>0):
            for j in range(len(indices_maxarea)):
                sliced = measurements.find_objects(lw == indices_maxarea[j])
                #sliced = measurements.find_objects(areaImg==maxarea)
                sliceX = sliced[0][1]
                sliceY = sliced[0][0]
                #print("len, sliced:", len(sliced))
                # Checking if the clusters span:
                if (sliceX.stop-sliceX.start)==Ly or (sliceY.stop-sliceY.start)==Lx:
                    '''
                    print("Spanning!:")
                    print("sliceX.start:", sliceX.start, ", sliceX.stop:", sliceX.stop)
                    print("sliceY.start:", sliceY.start, ", sliceY.stop:", sliceY.stop)
                    '''
                    P_av += values_maxarea[j]
                    #P_av += maxarea
                    break
                    '''
                    plotxlim=im3.axes.get_xlim()
                    plotylim=im3.axes.get_ylim()
                    plot([sliceX.start, sliceX.start, sliceX.stop, sliceX.stop, sliceX.start],
                    [sliceY.start, sliceY.stop, sliceY.stop, sliceY.start, sliceY.start],
                     color="red")
                    xlim(plotxlim)
                    ylim(plotylim)
                    show()
                    '''

    P_av /= (N_real*Lx*Ly)
    return P_av        

def findPi_av(Lx,Ly,p,N_real):
    Pi_av = 0
    p = float(p)
    for i in range(N_real):
        r = rand(Lx,Ly)
        z = r<p
        lw, num = measurements.label(z)
        
        b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
        shuffle(b) # shuffle this array
        shuffledLw = b[lw] # replace all values with values from b
        
        # Calculate areas
        area = measurements.sum(z, lw, index=arange(lw.max() + 1))
        areaImg = area[lw]
        # Bounding box
        # The cluster with the largest area is not always the one with the largest bounding box
        maxarea = areaImg.max()
        #indices_maxarea = indices_max_all(area)      # Finds the index of the cluster(s) that have the largest area
        indices_maxarea, value_maxarea = indices_Nlargest(area,10)   # Finds the index of the clusters that have areas corresponding to the N largest occurances

        # Just assuming that the largest cluster is the spanning cluster: (Seems fair from comparing with the other method)
        if(len(indices_maxarea)>0):
            for j in range(len(indices_maxarea)):
                sliced = measurements.find_objects(lw == indices_maxarea[j])
                #print("len, sliced:", len(sliced))
                if(len(sliced)>0):
                    sliceX = sliced[0][1]
                    sliceY = sliced[0][0]
                    #print("sliceX.start:", sliceX.start, ", sliceX.stop:", sliceX.stop)
                    #print("sliceY.start:", sliceY.start, ", sliceY.stop:", sliceY.stop)
                    # Checking if the cluster spans the system:
                    if (sliceX.stop-sliceX.start)==Lx or (sliceY.stop-sliceY.start)==Ly:
                        Pi_av += 1.0
                        break

    Pi_av /= N_real
    return Pi_av    

def findM_av(Lx,Ly,p,N_real):
    M_av = 0
    p = float(p)
    for i in range(N_real):
        r = rand(Lx,Ly)
        z = r<p
        #print(z)
        lw, num = measurements.label(z) # Labels cl
        
        #print(lw)
        b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
        shuffle(b) # shuffle this array
        shuffledLw = b[lw] # replace all values with values from b
        
        # Calculate areas
        area = measurements.sum(z, lw, index=arange(lw.max() + 1)) # Returns an array with the sum of all elements of label 1, 2, 3, etc., i.e. the areas, in the order specified by index
        areaImg = area[lw]
        
        '''
        figure()
        im3 = imshow(areaImg, origin='lower', interpolation='nearest')
        title("Clusters by area")
        '''
        
        #print(areaImg)
        #print(area)
        #print("area printed")
        #print("Type, area:", type(area))
        
        # Bounding box
        # The cluster with the largest area is not always the one with the largest bounding box
        maxarea = areaImg.max()
        #indices_maxarea = indices_max_all(area)      # Finds the index of the cluster(s) that have the largest area
        indices_maxarea, values_maxarea = indices_Nlargest(area,10)   # Finds the index of the clusters that have areas corresponding to the N largest occurances

        # Just assuming that the largest cluster is the spanning cluster: (Seems fair from comparing with the other method)
        counter = 0
        M_this  = 0
        if(len(indices_maxarea)>0):
            for j in range(len(indices_maxarea)):
                sliced = measurements.find_objects(lw == indices_maxarea[j]) # Find the object that has the same index as one of the largest areas
                #sliced = measurements.find_objects(areaImg==maxarea)
                sliceX = sliced[0][1]
                sliceY = sliced[0][0]
                #print("len, sliced:", len(sliced))
                if (sliceX.stop-sliceX.start)==Ly or (sliceY.stop-sliceY.start)==Lx:
                    '''
                    print("Spanning!:")
                    print("sliceX.start:", sliceX.start, ", sliceX.stop:", sliceX.stop)
                    print("sliceY.start:", sliceY.start, ", sliceY.stop:", sliceY.stop)
                    '''
                    M_this += values_maxarea[j]
                    counter+=1
                    
                    '''
                    plotxlim=im3.axes.get_xlim()
                    plotylim=im3.axes.get_ylim()
                    plot([sliceX.start, sliceX.start, sliceX.stop, sliceX.stop, sliceX.start],
                    [sliceY.start, sliceY.stop, sliceY.stop, sliceY.start, sliceY.start],
                     color="red")
                    xlim(plotxlim)
                    ylim(plotylim)
                    show()
                    '''
        if(counter!=0):
            M_av += M_this/counter
    M_av /= (N_real)
    return M_av  

def get_nps(pmin, pmax, N_ps, N_realizations, Lx, Ly, maxexp, base):
    if N_ps>1:
        nps = zeros((N_ps,maxexp+1))
        ps = linspace(pmin, pmax, N_ps)
    else:
        nps = zeros((N_ps,maxexp+1))
        ps = array([pmin])
    for i in range(N_realizations):                        # Looping over different realizations
        r = rand(Lx,Ly)
        for j in range(N_ps):                              # Looping over different p's
            z = r<ps[j]
            #print(z)
            lw, num = measurements.label(z)
        
            #print(lw)
            b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
            shuffle(b) # shuffle this array
            shuffledLw = b[lw] # replace all values with values from b
            
            # Calculate areas
            area = measurements.sum(z, lw, index=arange(lw.max() + 1))
            areaImg = area[lw]
            
            '''
            print("areaImg:",areaImg)
            print("area:",area)
            '''
            # Bounding box
            # The cluster with the largest area is not always the one with the largest bounding box
            maxarea = areaImg.max()
            #indices_maxarea = indices_max_all(area)      # Finds the index of the cluster(s) that have the largest area
            #print("area before finding indices:", area)
            indices_maxarea, value_maxarea = indices_Nlargest(area,10)   # Finds the index of the clusters that have areas corresponding to the N largest occurances
            area_indices_maxarea = copy(indices_maxarea)
            area_value_maxarea   = copy(value_maxarea)
            '''
            print(type(area))
            print("Area before removing anything:", area)
            print(areaImg)
            print(j)
            print(ps[j])
            
            print("indices_maxarea:", indices_maxarea)
            print("value_maxarea:", value_maxarea)
            '''
            
            k_area = 0 # index for area 
    
            # Just assuming that the largest cluster is the spanning cluster: (Seems fair from comparing with the other method)
            if(len(indices_maxarea)>0):
                for k in range(len(indices_maxarea)):
                    #print("k =",k)
                    #print("indices_maxarea[k]=",indices_maxarea[k])
                    sliced = measurements.find_objects(lw == indices_maxarea[k]) # Matching to elements in lw
                    #print("len, sliced:", len(sliced))
                    #print("sliced:", sliced)
                    if(value_maxarea[k]<Lx or value_maxarea[k]<Ly):  # Ignore clusters that cannot possibly span the system.
                        break
                    if(len(sliced)>0):
                        #print("In if-test, sliced=",sliced)
                        sliceX = sliced[0][1]
                        sliceY = sliced[0][0]
                        #print("sliceX.start:", sliceX.start, ", sliceX.stop:", sliceX.stop)
                        #print("sliceY.start:", sliceY.start, ", sliceY.stop:", sliceY.stop)
                        if (sliceX.stop-sliceX.start)==Lx or (sliceY.stop-sliceY.start)==Ly:
                            '''
                            print("Have spanning cluster")
                            print("sliceX.start:", sliceX.start, ", sliceX.stop:", sliceX.stop)
                            print("sliceY.start:", sliceY.start, ", sliceY.stop:", sliceY.stop)
                            print("area:", area)
                            print("len, area:", len(area))
                            print("area_indices_maxarea[k_area]:",area_indices_maxarea[k_area])
                            print("area[int(area_indices_maxarea[k_area])]=",area[int(area_indices_maxarea[k_area])])
                            print("area[]=",area[int(area_indices_maxarea[k_area])])
                            '''
                            area = delete(area, area_indices_maxarea[k_area])
                            area_indices_maxarea, area_value_maxarea = indices_Nlargest_delshift(indices_maxarea, value_maxarea,k) # Delete spanning cluster from area and shift indices
                            k_area -= 1
                    k_area+=1
            # Take the histogram
            maxarea_after = max(area)
            bc, histval = logarithmic_binning_setbins(maxexp,base,area,1)
            nps[j][:] = nps[j][:]+histval # Is this the format?
            #nps[j][:] = nps[j][:]+histval # Is this the format?
    # Take average of nps.
    nps /= (N_realizations*Lx*Ly)
    # Return nps. And maybe something else?
    return nps, bc, ps


def logarithmic_binning_setbins(largestexp, base, data, plot):
    #   Returning a histogram with exponentially increasing bin size
    #   Binning from zero and up to a set, maximal value
    #   Fed in by the user through the number of bins N
    #   and the base.
    #   Use the results from this with bar to plot

    # Finding the rightmost edge of each bin
    N = largestexp+1
    binedges = zeros(N)
    for i in range(N):
        binedges[i] = base**i # Could probably do this in vector form...
    
    # Finding the bin edges
    bin_length = zeros(N)
    bin_length[0] = 1
    for i in range(1,N):
        bin_length[i] = binedges[i]-binedges[i-1]

    bincenters    = binedges-bin_length/2.0 # To plot, maybe?
    binedges_left = binedges-bin_length
    loghist       = zeros(N)
    # Placing the data into bins
    for i in range(len(data)):
        thisdata = data[i]
        for j in range(N):
            if thisdata<=binedges[j] and thisdata>binedges_left[j]:
                loghist[j] = loghist[j] + 1
    for i in range(N):
        loghist[i] /= bin_length[i]
    #loghist = loghist/bin_length
    
    if plot==0:
        figure()
        #bar(bincenters,loghist) %draws the bars at the locations specified by x.
        loglog(bincenters,loghist)
        title('Logarithmic plotting')
        xlabel('Bin')
        ylabel('Value')
        show()
    return bincenters, loghist

### For MD data handling: ###

def is_the_given_matrix_percolating(vmat):
    isitpercolating = 0 # 0 if it isn't, 1 if it is
    isitpercolating_x = 0
    isitpercolating_y = 0
    isitpercolating_z = 0
    
    # Get dimensions
    dims      = np.shape(vmat) 
    Lx        = dims[0]
    Ly        = dims[1]
    Lz        = dims[2]
    print('Lx:', Lx, 'Ly:', Ly, 'Lz:', Lz)
    
    lw, num = measurements.label(vmat)
    
    print('num:', num)
    
    b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
    shuffle(b) # shuffle this array
    shuffledLw = b[lw] # replace all values with values from b
    
    # Calculate areas
    area = measurements.sum(vmat, lw, index=arange(lw.max() + 1))
    areaImg = area[lw]
    
    #print('Number of clusters:', len(area))
    # Bounding box
    # The cluster with the largest area is not always the one with the largest bounding box
    maxarea = areaImg.max()
    #indices_maxarea = indices_max_all(area)      # Finds the index of the cluster(s) that have the largest area
    indices_maxarea, value_maxarea = indices_Nlargest(area,10)   # Finds the index of the clusters that have areas corresponding to the N largest occurances
    
    # Just assuming that the largest cluster is the spanning cluster: (Seems fair from comparing with the other method)
    #print('NEW MATRIX!!')
    if(len(indices_maxarea)>0 and num!=0):
        for j in range(len(indices_maxarea)):
            sliced = measurements.find_objects(lw == indices_maxarea[j])
            #print("len, sliced:", len(sliced))
            if(len(sliced)>0):
                if num==1:
                    print('num =', num ,'; sliced:',sliced)
                sliceX = sliced[0][0]  # I have no idea if this works in 3d...
                sliceY = sliced[0][1]
                sliceZ = sliced[0][2]
                #print("sliceX.start:", sliceX.start, ", sliceX.stop:", sliceX.stop)
                #print("sliceY.start:", sliceY.start, ", sliceY.stop:", sliceY.stop)
                # Checking if the cluster spans the system:
                if (sliceX.stop-sliceX.start)==Lx:
                    isitpercolating_x = 1.0
                if (sliceY.stop-sliceY.start)==Ly:
                    isitpercolating_y = 1.0
                if (sliceZ.stop-sliceZ.start)==Lz:
                    isitpercolating_z = 1.0
                if isitpercolating_x==1.0 and isitpercolating_y==1.0 and isitpercolating_z==1.0:
                    isitpercolating = 1.0
                    break
    return isitpercolating, isitpercolating_x, isitpercolating_y, isitpercolating_z


if __name__=="__main__":
    Lx = 100
    Ly = 5#Lx
    p = 0.65
    N_real = 10#00
    
    P_av = findP_av(Lx,Ly,p,N_real)
    print("p=",p)
    print(P_av)
    
    largestexp = 10
    base = 2
    plot = 0
    data = numpy.random.rand(1,100000)
    data = data[0]
    for i in range(len(data)):
        data[i] = 1/(data[i]**2)
    #print(data)
    #bc, loghist =logarithmic_binning_setbins(largestexp, base, data, plot)
    #print(loghist)
    
    pmin =0.2; pmax=0.9; N_ps=10; N_realizations=10; Lx=10; Ly=10; maxexp=10; base=2;
    nps = get_nps(pmin, pmax, N_ps, N_realizations, Lx, Ly, maxexp, base)
    
