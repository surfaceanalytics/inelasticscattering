# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:46:57 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt

from data_converter import DataConverter

class DataProcessor(DataConverter):
    def __init__(self):
        super().__init__()
           
    def filterRegion(self, region):
        """ This method performs a query on the data set, and returns the
        only the data for which the 'group.region.name' key is equal to the 
        inputted 'region' variable.
        """
        data_list = []    
        for _group in self.data['groups']:
            for _region in _group['regions']:
                if _region['name'] == region:
                    for n, cycle in enumerate(_region['cycles']):
                        data_list += [cycle['data']]
        return data_list
            
    def listRegions(self):
        """ This method returns a list of all the regions names in the dataset.
        """
        regions_list = []
        for _group in self.data['groups']:
            for _region in _group['regions']:
                regions_list += [_region['name']]
        return regions_list
                
    def getAllRegions(self):
        """ This method returns a list dictionary of dictionaries, where the
        top-most key is the region name (i.e. 'O1s', 'C1s', etc.), and the
        values are a list of all spectra having this region name.
        """
        all_regions = {}
        for i in self.listRegions():
            all_regions[i] = self.filterRegion(i)
        return all_regions
            
    def normalize(self):
        """ This method divides the main data channel (i.e. 'y0') by any 
        additional data channels (if they exist).
        """
        normalized = {}
        for idx, region in self.getAllRegions().items():
            normalized[idx] = {}
            for cycle, data in region.items():
                normalized[idx][cycle] = {}
                normalized[idx][cycle]['x'] = data['x']                
                normalized[idx][cycle]['y0'] = data['y0']
                for channel in range(self.nr_channels-1):
                    normalized[idx][cycle]['ny'+str(channel+1)] = list(np.divide(data['y0'],data['y'+str(channel+1)]))
        return normalized
        
#%%
data = DataProcessor() 
filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\XPS data analysis manuscript with Neal\C1s scans.xy'       
data.load(filepath)
data.write('test', out_format = 'Vamas')

# Get only the spectra from one region
# Returns the data in a dictionary
c1s = data.filterRegion('C1s')

# Get only the spectra from all regions
# Returns the data in a dictionary
all_regions = data.getAllRegions()

normalized = data.normalize()

#%%

# this gets just the y values from the data and puts it in a matrix
data_collection = np.matrix([data.data[0]['regions'][1]['cycles'][i]['data']['y0']
    for i in data.data[0]['regions'][1]['cycles']])

# This sums together all previsous spectra for 
sums = np.matrix(np.array([np.sum(data_collection[:i], axis=0) 
    for i in range(len(data_collection))]))[1:]

# This sums in a different way
sum_of_counts = np.matrix(np.zeros(np.shape(data_collection)[1]))
persist = np.matrix(np.zeros(np.shape(data_collection)[1]))
for i in data_collection:
    persist = np.sum([i,persist], axis=0)
    sum_of_counts = np.append(sum_of_counts,persist, axis=0)

sum_of_counts = sum_of_counts[1:,:]

# this divides the sum of counts by the number of scans
counts_per_sec = np.array([sum_of_counts[i] / (i+1) 
    for i in range(len(sum_of_counts[:,0]))])[:,0,:]

# this shows the first and last scan
plt.plot(counts_per_sec[0,:].T)
plt.plot(counts_per_sec[-1,:].T)
plt.show()


# this gets the areas under each curve
areas = []
for i in counts_per_sec:
    areas+=[np.average(i)]

plt.plot(areas[:800])


# This subtracts the last scan (taken as an approximation to the true
# signal) from each of the scans.
difference = []
for i in counts_per_sec:
    difference += [np.subtract(i, counts_per_sec[-1,:])]

difference = np.array(difference)

plt.plot(difference[0,:])
plt.plot(difference[-200,:])
plt.show()


# This sums together the intensity of the noise
noise = []
for i in difference:
    positive = np.absolute(i) 
    noise += [np.sum(positive)]
    
plt.plot(noise)  


# This calculates the signal-noise ratio

# the vriable 'signal' is taken to be the true signal intensity
signal = areas[-1] - np.min(counts_per_sec[-1]) 
sig_to_noise = [signal / i for i in noise]

plt.plot(sig_to_noise[:-500])
 

#%%

# Put the averaged data back into the data object
data_to_write = sum_of_counts

for i in range(len(data_to_write[:,0])):
    data.data[0]['regions'][1]['cycles'][i]['data']['y0'] = [d for d in data_to_write[i,:].T]

# export to vamas
data.writeVamas('loops_sum_of_counts')
