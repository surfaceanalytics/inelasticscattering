# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 17:29:54 2020

@author: Mark
"""
from vamas import VamasHeader, Block

class VamasParser():
    def __init__(self):
        self.header = VamasHeader()
        self.blocks = []
        self.common_header_attr = ['formatID', 'instituteID', 
                                  'instrumentModelID', 'operatorID',
                                  'experimentID', 'noCommentLines']
        
        self.norm_header_attr = ['scanMode','unknown1','unknown2', 
                                 'expVarLabel','expVarUnit','unknown3',
                                 'unknown4','unknown5','unknown6','noBlocks']
        
        self.map_header_attr = ['scanMode','unknown1','nr_positions',
                                'nr_x_coords', 'nr_y_coords','unknown2', 
                                 'expVarLabel','expVarUnit','unknown3',
                                 'unknown4','unknown5','unknown6','noBlocks']
        
        self.norm_block_attr = ['blockID', 'sampleID', 'year', 'month', 'day',
                                'hour', 'minute', 'second', 
                                'noHrsInAdvanceOfGMT', 'noCommentLines', 
                                'commentLines', 'technique', 'expVarValue', 
                                'sourceLabel', 'sourceEnergy', 'unknown1', 
                                'unknown2', 'unknown3', 'sourceAnalyzerAngle', 
                                'unknown4', 'analyzerMode', 'resolution', 
                                'magnification', 'workFunction', 'targetBias', 
                                'analyzerWidthX', 'analyzerWidthY', 
                                'analyzerTakeOffPolarAngle', 'analyzerAzimuth', 
                                'speciesLabel', 'transitionLabel', 
                                'particleCharge', 'abscissaLabel', 
                                'abscissaUnits', 'abscissaStart', 
                                'abscissaStep', 'noVariables', 'variableLabel1', 
                                'variableUnits1', 'variableLabel2', 
                                'variableUnits2', 'signalMode', 'dwellTime', 
                                'noScans', 'timeCorrection', 'sampleAngleTilt', 
                                'sampleTiltAzimuth', 'sampleRotation', 
                                'noAdditionalParams', 'paramLabel1', 'paramUnit1', 
                                'paramValue1', 'paramLabel2', 'paramUnit2', 
                                'paramValue2', 'numOrdValues', 'minOrdValue1', 
                                'maxOrdValue1', 'minOrdValue2', 'maxOrdValue2', 
                                'dataString']
        
        self.map_block_attr = ['blockID', 'sampleID', 'year', 'month', 'day',
                                'hour', 'minute', 'second', 
                                'noHrsInAdvanceOfGMT', 'noCommentLines', 
                                'commentLines', 'technique', 'x_coord', 
                                'y_coord', 'expVarValue', 
                                'sourceLabel', 'sourceEnergy', 'unknown1', 
                                'unknown2', 'unknown3', 'fov_x', 'fovy',
                                'sourceAnalyzerAngle', 'unknown4', 
                                'analyzerMode', 'resolution', 'magnification', 
                                'workFunction', 'targetBias', 'analyzerWidthX', 
                                'analyzerWidthY', 'analyzerTakeOffPolarAngle', 
                                'analyzerAzimuth', 'speciesLabel', 'transitionLabel', 
                                'particleCharge', 'abscissaLabel', 
                                'abscissaUnits', 'abscissaStart', 
                                'abscissaStep', 'noVariables', 'variableLabel1', 
                                'variableUnits1', 'variableLabel2', 
                                'variableUnits2', 'signalMode', 'dwellTime', 
                                'noScans', 'timeCorrection', 'sampleAngleTilt', 
                                'sampleTiltAzimuth', 'sampleRotation', 
                                'noAdditionalParams', 'paramLabel1', 'paramUnit1', 
                                'paramValue1', 'paramLabel2', 'paramUnit2', 
                                'paramValue2', 'numOrdValues', 'minOrdValue1', 
                                'maxOrdValue1', 'minOrdValue2', 'maxOrdValue2', 
                                'dataString']

    def parseFile(self, filepath):
        """ The openFile method parses the .xy file into a list of dictionaries 
        under the attribute 'self.data'.
        Each dictionary is a grouping of related attributes.
        These are later put into a heirarchical nested dictionary that 
        represents the native data structure of the export, and is well 
        represented by JSON.
        """
        self._readLines(filepath)
        self._parseHeader()
        self._parseBlocks()
        return self._buildDict()

    def _readLines(self, filepath):
        self.data = []
        self.filepath = filepath
        with open(filepath) as fp:
            for line in fp:
               self.data += [line]
        
    def _parseHeader(self):
        for attr in self.common_header_attr:
            setattr(self.header, attr, self.data.pop(0).strip())
        n = int(self.header.noCommentLines)
        comments = ''
        for l in range(n):
            comments += self.data.pop(0)
        self.header.commentLines = comments
        self.header.expMode = self.data.pop(0).strip()
        if self.header.expMode == 'NORM':
            for attr in self.norm_header_attr:
                setattr(self.header, attr, self.data.pop(0).strip())
        elif self.header.expMode == 'MAP':
            for attr in self.map_header_attr:
                setattr(self.header, attr, self.data.pop(0).strip())
        
    def _parseBlocks(self):
        for b in range(int(self.header.noBlocks)):
            self._parseOneBlock()
    
    def _parseOneBlock(self):
        if self.header.expMode == 'NORM':
            self.blocks += [self._parseNORMBlock()]
        elif self.header.expMode == 'MAP':
            self.blocks += [self._parseMAPBlock()]
            
    def _parseNORMBlock(self):
        block = Block()
        block.blockID = self.data.pop(0).strip()
        block.sampleID = self.data.pop(0).strip()
        block.year = int(self.data.pop(0).strip())
        block.month = int(self.data.pop(0).strip())
        block.day = int(self.data.pop(0).strip())
        block.hour = int(self.data.pop(0).strip())
        block.minute = int(self.data.pop(0).strip())
        block.second = int(self.data.pop(0).strip())
        block.noHrsInAdvanceOfGMT = int(self.data.pop(0).strip())
        block.noCommentLines = int(self.data.pop(0).strip())
        for n in range(block.noCommentLines):
            block.commentLines += self.data.pop(0)       
        block.technique = self.data.pop(0).strip()
        block.expVarValue = self.data.pop(0).strip()
        block.sourceLabel = self.data.pop(0).strip()
        block.sourceEnergy = float(self.data.pop(0).strip())
        block.unknown1 = self.data.pop(0).strip()
        block.unknown2 = self.data.pop(0).strip()
        block.unknown3 = self.data.pop(0).strip()
        block.sourceAnalyzerAngle = self.data.pop(0).strip()
        block.unknown4 = self.data.pop(0).strip()
        block.analyzerMode = self.data.pop(0).strip()
        block.resolution = float(self.data.pop(0).strip())
        block.magnification = self.data.pop(0).strip()
        block.workFunction = float(self.data.pop(0).strip())
        block.targetBias = float(self.data.pop(0).strip())
        block.analyzerWidthX = self.data.pop(0).strip()
        block.analyzerWidthY = self.data.pop(0).strip()
        block.analyzerTakeOffPolarAngle = self.data.pop(0).strip()
        block.analyzerAzimuth = self.data.pop(0).strip()
        block.speciesLabel = self.data.pop(0).strip()
        block.transitionLabel = self.data.pop(0).strip()
        block.particleCharge = self.data.pop(0).strip()
        block.abscissaLabel = self.data.pop(0).strip()
        block.abscissaUnits = self.data.pop(0).strip()
        block.abscissaStart = float(self.data.pop(0).strip())
        block.abscissaStep = float(self.data.pop(0).strip())
        block.noVariables = int(self.data.pop(0).strip())
        for p in range(block.noVariables):
            name = 'variableLabel' + str(p+1)
            setattr(block, name, self.data.pop(0).strip())
            name = 'variableUnits' + str(p+1)
            setattr(block, name, self.data.pop(0).strip())
        block.signalMode = self.data.pop(0).strip()
        block.dwellTime = float(self.data.pop(0).strip())
        block.noScans = int(self.data.pop(0).strip())
        block.timeCorrection = self.data.pop(0).strip()
        block.sampleAngleTilt = float(self.data.pop(0).strip())
        block.sampleTiltAzimuth = float(self.data.pop(0).strip())
        block.sampleRotation = float(self.data.pop(0).strip())
        block.noAdditionalParams = int(self.data.pop(0).strip())
        for p in range(block.noAdditionalParams):
            name = 'paramLabel' + str(p+1)
            setattr(block, name, self.data.pop(0))
            name = 'paramUnit' + str(p+1)
            setattr(block, name, self.data.pop(0))
            name = 'paramValue' + str(p+1)
            setattr(block, name, self.data.pop(0))
        block.numOrdValues = int(self.data.pop(0).strip())
        for p in range(block.noVariables):
            name = 'minOrdValue' + str(p+1)
            setattr(block, name, float(self.data.pop(0).strip()))
            name = 'maxOrdValue' + str(p+1)
            setattr(block, name, float(self.data.pop(0).strip()))
            
        self._addDataValues(block)

        return block
    
    def _parseMAPBlock(self):
        block = Block()
        block.blockID = self.data.pop(0).strip()
        block.sampleID = self.data.pop(0).strip()
        block.year = int(self.data.pop(0).strip())
        block.month = int(self.data.pop(0).strip())
        block.day = int(self.data.pop(0).strip())
        block.hour = int(self.data.pop(0).strip())
        block.minute = int(self.data.pop(0).strip())
        block.second = int(self.data.pop(0).strip())
        block.noHrsInAdvanceOfGMT = int(self.data.pop(0).strip())
        block.noCommentLines = int(self.data.pop(0).strip())
        for n in range(block.noCommentLines):
            block.commentLines += self.data.pop(0)       
        block.technique = self.data.pop(0).strip()
        block.x_coord = self.data.pop(0).strip()
        block.y_coord = self.data.pop(0).strip()
        block.expVarValue = self.data.pop(0).strip()
        block.sourceLabel = self.data.pop(0).strip()
        block.sourceEnergy = float(self.data.pop(0).strip())
        block.unknown1 = self.data.pop(0).strip()
        block.unknown2 = self.data.pop(0).strip()
        block.unknown3 = self.data.pop(0).strip()
        block.fov_x = self.data.pop(0).strip()
        block.fov_y = self.data.pop(0).strip()
        block.sourceAnalyzerAngle = self.data.pop(0).strip()
        block.unknown4 = self.data.pop(0).strip()
        block.analyzerMode = self.data.pop(0).strip()
        block.resolution = float(self.data.pop(0).strip())
        block.magnification = self.data.pop(0).strip()
        block.workFunction = float(self.data.pop(0).strip())
        block.targetBias = float(self.data.pop(0).strip())
        block.analyzerWidthX = self.data.pop(0).strip()
        block.analyzerWidthY = self.data.pop(0).strip()
        block.analyzerTakeOffPolarAngle = self.data.pop(0).strip()
        block.analyzerAzimuth = self.data.pop(0).strip()
        block.speciesLabel = self.data.pop(0).strip()
        block.transitionLabel = self.data.pop(0).strip()
        block.particleCharge = self.data.pop(0).strip()
        block.abscissaLabel = self.data.pop(0).strip()
        block.abscissaUnits = self.data.pop(0).strip()
        block.abscissaStart = float(self.data.pop(0).strip())
        block.abscissaStep = float(self.data.pop(0).strip())
        block.noVariables = int(self.data.pop(0).strip())
        for p in range(block.noVariables):
            name = 'variableLabel' + str(p+1)
            setattr(block, name, self.data.pop(0).strip())
            name = 'variableUnits' + str(p+1)
            setattr(block, name, self.data.pop(0).strip())
        block.signalMode = self.data.pop(0).strip()
        block.dwellTime = float(self.data.pop(0).strip())
        block.noScans = int(self.data.pop(0).strip())
        block.timeCorrection = self.data.pop(0).strip()
        block.sampleAngleTilt = float(self.data.pop(0).strip())
        block.sampleTiltAzimuth = float(self.data.pop(0).strip())
        block.sampleRotation = float(self.data.pop(0).strip())
        block.noAdditionalParams = int(self.data.pop(0).strip())
        for p in range(block.noAdditionalParams):
            name = 'paramLabel' + str(p+1)
            setattr(block, name, self.data.pop(0))
            name = 'paramUnit' + str(p+1)
            setattr(block, name, self.data.pop(0))
            name = 'paramValue' + str(p+1)
            setattr(block, name, self.data.pop(0))
        block.numOrdValues = int(self.data.pop(0).strip())
        for p in range(block.noVariables):
            name = 'minOrdValue' + str(p+1)
            setattr(block, name, float(self.data.pop(0).strip()))
            name = 'maxOrdValue' + str(p+1)
            setattr(block, name, float(self.data.pop(0).strip()))
            
        self._addDataValues(block)
        
        return block
    
    def _addDataValues(self, block):
        data_dict = {}
        start = float(block.abscissaStart)
        step = float(block.abscissaStep)
        num = int(block.numOrdValues / block.noVariables)
        x = [round(start + i*step,2) for i in range(num)]
        
        if block.abscissaLabel == 'binding energy':
            x.reverse()
        
        setattr(block, 'x', x) 
        
        for v in range(block.noVariables):
            name = 'y' + str(v)
            data_dict[name] = [] 
            
        for r in range(int(block.numOrdValues / block.noVariables)):
            for v in range(block.noVariables):
                name = 'y' + str(v)
                data_dict[name] += [float(self.data.pop(0).strip())]
                
        for v in range(block.noVariables):
            name = 'y' + str(v)
            setattr(block, name, data_dict[name]) 

    def _buildDict(self):
        """ This method constructs a list of dictionaries where eache 
        dictionary contains all the data and metadata of a spectrum.
        vamas.sampleID -> group['name']
        vamas.
        """
        
        group_id = -1
        temp_group_name = ''
        spectra = []
        
        for idx, b in enumerate(self.blocks):
            group_name = b.sampleID
            
            ''' This set of conditions detects if the group name has changed.
            If it has, then it increments the group_idx.
            '''
            if group_name != temp_group_name:
                temp_group_name = group_name
                group_id += 1
           
            spectrum_type = str(b.speciesLabel + b.transitionLabel)
            spectrum_id = idx   
            
            settings = {
                    'analysis_method':b.technique,
                    'dwell_time':b.dwellTime,
                    'workfunction':b.workFunction,
                    'excitation_energy':b.sourceEnergy,
                    'pass_energy':b.resolution,
                    'scan_mode':b.analyzerMode,
                    'source_label':b.sourceLabel,
                    'nr_values':int(b.numOrdValues / b.noVariables),
                    'x_units': b.abscissaLabel,
                    'y_units': b.variableLabel1
                    }
            
            date = (str(b.year) + '-' + str(b.month) + '-' + str(b.day) 
                    + ' ' + str(b.hour) + ':' + str(b.minute) + ':' 
                    + str(b.second))
                
            data = {'x':b.x}
            for n in range(int(b.noVariables)):
                key = 'y'+str(n)
                data[key] = getattr(b, key)

            spec_dict = {'date':date,
                         'group_name': group_name, 
                         'group_id':group_id,
                         'spectrum_type':spectrum_type,
                         'spectrum_id':spectrum_id,
                         'scans':b.noScans, 
                         'settings':settings, 
                         'data':data}
                
            spectra += [spec_dict]

        self.data_dict = spectra
        return self.data_dict

#%%
if __name__ == '__main__':
    filepath = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Data Science\xps_data_conversion_tools\EX337 - test.vms'
    v = VamasParser()
    v.parseFile(filepath)
    