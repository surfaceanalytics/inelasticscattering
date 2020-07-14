# -*- coding: utf-8 -*-
"""
Created on Tue Jun 2 22:11:04 2020

@author: Mark
"""

class Dataset():
    """ The 'Dataset' object is simply a container for the data, to be used for 
    converting spectra, with its experimental metadata, into a different 
    format, such as vamas.
    """
    def __init__(self):
        self.num_spectra = 0
        self.blocks = []
        self.sourceAnalyzerAngle = '56.5'
        self.vamas_header = VamasHeader()
        
    def data_to_blocks(self, data):
        """ This converts the 'data' which should be a dictionary of 
        dictionaries into a 'Block'. A 'Block' object is used in Vamas.
        It is essentially the raw data plus all experimental metadata for a
        spectrum. There is one block per spectrum.
        The attributes of the data are set to Block attributes.
        The blocks are then stored in a list of blocks called self.block.
        """
        idx=0
        for d in data:
            block = Block()
            block.blockID = d['blockID']
            block.sampleID = d['sampleID']
            block.month = d['month']
            block.day = d['day']
            block.year = d['year']
            block.hour = d['hour']
            block.minute = d['minute']
            block.second = d['second'] 
            block.noHrsInAdvanceOfGMT = '0'
            block.noCommentLines = 5
            block.commentLines = 'Casa Info Follows\n0\n0\n0\n0'
            block.technique = "XPS"
            block.expVarValue = 0
            block.sourceLabel = 'Al'
            block.sourceEnergy = 1486.7
            block.sourceAnalyzerAngle = 0
            block.analyzerMode = 'FAT'
            block.resolution = 0
            block.workFunction = 0
            block.speciesLabel = 'species'
            block.transitionLabel = 'transition'
            block.abscissaStart = d['x'][0]
            block.abscissaStep = round(abs(d['x'][0]-d['x'][1]),6)
            block.dwellTime = 0.1
            block.noScans = 1
            block.numOrdValues = int(len(d['y'])*2)
            block.minOrdValue1 = round(min(d['y']),6)
            block.maxOrdValue1 = round(max(d['y']),6)
            block.minOrdValue2 = 1
            block.maxOrdValue2 = 1
            for i in d['y']:
                ''' This for loop appends the ordinate data, i.e. the y values,
                to the data stream. It separates the vales by a line break,
                and interlaaces a 1 in between each ordinate value. The 1 is
                used because, in common use, there are two "channels" in the 
                data stream: one for the spectrum data, and one for the
                transmission function. Since the transmission function is 
                not in our case useful, we leave it as a constant 1.
                '''
                block.dataString += str(round(i,6)) + '\n1\n'
            block.dataString = block.dataString[:-1]
            self.blocks += [block]
            idx += 1
            
        ''' The vamas header needs to state the number of blocks in the data
        file. So here we update the .noBlocks attribute of the vamas_header.
        '''    
        self.num_spectra = len(self.blocks)
        self.vamas_header.noBlocks = self.num_spectra
                           
    def writeVamas(self, filename):
        """ This writes the vamas header and blocks to an ASCII encoded file
        with the extension .vms, i.e. a vamas file.
        It first iterates through the vamas_header attributes, writes them to
        ASCII, and adds a line break after every attribute.
        Then it goes through the list of blocks, and for each block,
        it iterates through the block's attributes, writing them to ASCII,
        and adding a line break at the end of each
        """
        with open(str(filename), 'w') as file:
            for item in self.vamas_header.__dict__:    
                file.writelines(str(self.vamas_header.__dict__[item]) + '\n')
            for block in self.blocks:
                for item in block.__dict__:
                    file.writelines(str(block.__dict__[item]) + '\n')
            file.writelines('end of experiment')
            file.close()
        
class VamasHeader():
    """ This object is to hold the header, which needs to be present at the
    top of a vamas file. The order in which these values are defines is 
    important!
    """
    def __init__(self):
        self.formatID = 'VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4'
        self.instituteID = 'Not Specified'
        self.instriumentModelID = 'Not Specified'
        self.operatorID = 'Not Specified'
        self.experimentID = 'Not Specified'
        self.noCommentLines = '2'
        self.commentLines = 'Casa Info Follows CasaXPS Version 2.3.22PR1.0\n0'
        self.expMode = 'NORM'
        self.scanMode = 'REGULAR'
        self.unknown1 = '0'
        self.unknown2 = '1'
        self.expVarLabel = 'Exp Variable'
        self.expVarUnit = 'd'
        self.unknown3 = '0'
        self.unknown4 = '0'
        self.unknown5 = '0'
        self.unknown6 = '0'
        self.noBlocks = '1'
        
class Block():
    """ A Block is used to store all the spectrum data and corresponding 
    experimental metadata for a spectrum. Generally there is one block 
    for each spectrum. 
    The order that the attributes are defined is important!
    """
    def __init__(self):
        self.blockID = ''
        self.sampleID = ''
        self.year = ''
        self.month = ''
        self.day = ''
        self.hour = ''
        self.minute = ''
        self.second = ''
        self.noHrsInAdvanceOfGMT = '0'
        self.noCommentLines = ''
        self.commentLines = '' #This list should contain one element for each line in the comment block
        self.technique = ''
        self.expVarValue = ''
        self.sourceLabel = ''
        self.sourceEnergy = ''
        self.unknown1 = '0'
        self.unknown2 = '0'
        self.unknown3 = '0'
        self.sourceAnalyzerAngle = ''
        self.unknown4 = '180'
        self.analyzerMode = ''
        self.resolution = ''
        self.magnification = '1'
        self.workFunction = ''
        self.targetBias = '0'
        self.analyzerWidthX = '0'
        self.analyzerWidthY = '0'
        self.analyzerTakeOffPolarAngle = '0'
        self.analyzerAzimuth = '0'
        self.speciesLabel = ''
        self.transitionLabel = ''
        self.particleCharge = '-1'
        self.abscissaLabel = 'kinetic energy'
        self.abscissaUnits = 'eV'
        self.abscissaStart = ''
        self.abscissaStep = ''
        self.noVariables = '2'
        self.variableLabel1 = 'counts'
        self.variableUnits1 = 'd'
        self.variableLabel2 = 'Transmission'
        self.variableUnits2 = 'd'
        self.signalMode = 'pulse counting'
        self.dwellTime = ''
        self.noScans = ''
        self.timeCorrection = '0'
        self.sampleAngleTilt = '0'
        self.sampleTiltAzimuth = '0'
        self.sampleRotation = '0'
        self.noAdditionalParams = '2'
        self.paramLabel1 = 'ESCAPE DEPTH TYPE'
        self.paramUnit1 = 'd'
        self.paramValue1 = '0'
        self.paramLabel2 = 'MFP Exponent'
        self.paramUnit2 = 'd'
        self.paramValue2 = '0'
        self.numOrdValues = ''
        self.minOrdValue1 = ''
        self.maxOrdValue1 = ''
        self.minOrdValue2 = ''
        self.maxOrdValue2 = ''
        self.dataString = ''
        



    