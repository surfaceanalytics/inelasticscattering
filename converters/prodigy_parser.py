# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:05:25 2020

@author: Mark
"""

from datetime import datetime

class ProdigyParser():
    """This class is used for parsing ASCII-encoded data exported from 
    SpecsLab Prodigy v 4.64.1-r88350.
    The native ASCII export file type is '.xy'.
    """
    def __init__(self, **kwargs): 
        self.normalize = 0
        if 'n_channels' in kwargs.keys():
            self.nr_channels = kwargs['n_channels']
        else:
            self.nr_channels = 0
            
        self.settings_map = {'Acquisition Date':'date',
                             'Analysis Method':'analysis_method',
                             'Analyzer Lens':'lens',
                             'Analyzer Slit': 'entrance_slit',
                             'Bias Voltage':'bias_voltage',
                             'Binding Energy':'binding_energy',
                             'Scan Mode': 'scan_mode',
                             'Values/Curve':'nr_values',
                             'Eff. Workfunction':'workfunction',
                             'Excitation Energy':'excitation_energy',
                             'Dwell Time':'dwell_time',
                             'Detector Voltage':'detector_voltage',
                             'Comment':'comments',
                             'Curves/Scan':'curves_per_scan',
                             'Pass Energy':'pass_energy',
                             'Source':'source_label'}
             
    def parseFile(self, filepath, commentprefix = '#', headerlines = 15):
        """ The openFile method parses the .xy file into a list of dictionaries 
        under the attribute 'self.data'.
        Each dictionary is a grouping of related attributes.
        These are later put into a heirarchical nested dictionary that 
        represents the native data structure of the export, and is well 
        represented by JSON.
        """
        self.data = {}
        self.prefix = commentprefix
        self.filepath = filepath
        self.headerlines = headerlines
        self.header = {}     
        with open(filepath) as fp:
            nr = 0
            block = 0
            last_number = False
            open_block = False
            channel_count = 0
            self.data[block] = {}
            for line in fp:
                # This first parses the header
                if nr < self.headerlines:
                    line = line.strip(self.prefix).strip()
                    if len(line) == 0:
                        pass
                    else:
                        self.header[line.split(':',1)[0].strip()] = line.split(':',1)[1].strip()
                # if the line starts with a comment prefix and that is the 
                # only thing in the line or if the line is completely empty
                # then close the block
                else:
                    # This intializes a new block
                    if (
                            ((line[0] == self.prefix) & 
                            (len(line.strip(self.prefix).strip()) == 0)) | 
                            (len(line.strip())==0) |
                            ((last_number == True) & (line[0] == self.prefix))
                            ):
                        open_block = False
                        last_number = False
                    # This if statement parses the metadata
                    # if the line starts with the prefix, and the prefix is not 
                    # the only thing in the line then check if a block is 
                    # already open. If not, then start a new block then add 
                    # data to the block
                    if((line[0] == self.prefix) & 
                            (len(line.strip(self.prefix).strip()) != 0)
                            ):
                        line = line.strip(self.prefix).strip()
                        if open_block == False:
                            open_block = True
                            block += 1
                            self.data[block]={}
                        key = line.split(':',1)[0]
                        val = line.split(':',1)[1].strip()
                        self.data[block][key] = val
                        # this is to count the number of channels in the data streams
                        # the word 'Cycle' indicates the start of a new set of data channels
                        if key == 'Cycle':
                            channel_count = 0
                            if (len(val.split(',')) > 1) and (':' in val):
                                comma_split = val.split(',')
                                self.data[block][key] = comma_split[0]
                                for s in comma_split[1:]:
                                    k = s.split(':')[0].strip()
                                    v = s.split(':')[1].strip()
                                    self.data[block][k] = v
                                
                    # if the line does not start with the prefix, and it is not empty
                    # then the line must start with a number
                    elif (line[0] != self.prefix) & (len(line.strip()) != 0):
                        # if no block is open yet, then open one
                        if open_block == False:
                            block += 1
                            channel_count +=1
                            # this is to indicate a new data stream. each data stream counts as a channel
                            if channel_count > self.nr_channels:
                                self.nr_channels = channel_count
                            self.data[block]={}
                            self.data[block]['data']=[]
                            open_block = True
                        self.data[block]['data'] += [[float(line.split(' ')[0]),float(line.split(' ')[-1])]]
                        last_number = True
                nr += 1
            self.data = [val for val in self.data.values() if len(val) !=0]
        return self._buildDict()         
    
    def _buildDict(self):
        """ This method constructs a nested dictionary that represent the native
        heirarchical data model in a form that is easily exported to JSON.
        """
        def rotate(l):
            l = l[1:]+l[:1]
            return l
        group_count = 0
        spectrum_count = 0
        channel = 1
        spectra = []
        data = {}
        for i in self.data:
            if 'Group' in i.keys():
                # add a new 'group' to the array of groups
                group_name = i['Group']
                group_id = group_count
                group_count += 1
            if 'Region' in i.keys():
                ''' Regions is a grouping of data, that is not used in Vamas.
                It is a form of normalization, but is not general purpose, so
                it is not used in the nested dictionary structure.'''
                spectrum_type = i['Region']
                settings = {key : i[key] for key in i.keys() if key != 'Region'}
                settings = self._replaceKeys(settings, self.settings_map)
                settings['y_units'] = self.header['Count Rate']
                settings['x_units'] = self.header['Energy Axis']
                if 'date' in settings.keys():
                     self._parseDatetime(settings['date'])
                
                n_scans = 1
                
            if 'Number of Scans' in i.keys():
                n_scans = int(i['Number of Scans'])
                
            if 'Acquisition Date' in i.keys():
                date = self._parseDatetime(i['Acquisition Date'])
                date = date.strftime('%Y-%m-%d %H:%M:%S')

            if 'data' in i.keys():
                spectrum_id = spectrum_count
                
                if channel <= self.nr_channels:
                    data['x'] = [j[0] for j in i['data']]
                    data['y'+str(channel-1)] = [j[1] for j in i['data']]
                    if channel == self.nr_channels:
                        new_spectrum = {'date': date,
                                        'group_name': group_name,
                                        'group_id':group_id,
                                        'spectrum_type': spectrum_type,
                                        'spectrum_id': spectrum_id,
                                        'scans': n_scans,
                                        'settings':settings.copy(),
                                        'data':data.copy()}
                        spectra += [new_spectrum.copy()]
                        channel = 1
                        spectrum_count += 1
                    else: 
                        channel += 1

        self.data_dict = spectra
        return self.data_dict
    
    def _parseDatetime(self, date):
        #difference = 0
        if date.find('UTC'):
            #difference = int(date[date.find('UTC')+3:])
            date = date[:date.find('UTC')].strip()
        else:
            date = date.strip()
            
        date_object = datetime.strptime(date, '%m/%d/%y %H:%M:%S')
    
        return date_object
    
    def _replaceKeys(self, dictionary, key_map):
        for key in key_map.keys():
            if key in dictionary.keys():
                dictionary[key_map[key]] = dictionary[key]
                dictionary.pop(key, None)
        return dictionary
        

            
