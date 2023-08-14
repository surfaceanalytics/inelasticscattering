# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:32:48 2020

@author: Mark
"""

from converters.prodigy_parser import ProdigyParser
from converters.vamas_parser import VamasParser
from converters.text_parser import TextParser
from converters.writers import JSONWriter, VamasWriter, ExcelWriter

class DataConverter():
    """ This class parses files of type:
            'ProdigyXY', 'Vamas'
        And writes files of type:
            'JSON', 'Vamas', 'Excel'
    """

    def __init__(self):
        """ All objects of the Dataset class should have the methods specificed
        in the attribute 'self.class_methods'
        """
        self._parser_methods = {'ProdigyXY': ProdigyParser,
                                'Vamas': VamasParser,
                                'Text':TextParser}
        self._write_methods = {'JSON':JSONWriter,
                               'Vamas':VamasWriter,
                               'Excel':ExcelWriter}
        self._extensions = {'xy': 'ProdigyXY',
                            'vms': 'Vamas',
                            'txt':'Text'}

    def load(self, filename, **kwargs):
        """ This method parses an input file an places it into a nested
        dictionary.
        Parameters
        ----------
        filename: STRING
            The location and name of the file you wish to parse.
        **kwargs:
            in_format: The file format of the loaded file.
        """
        if 'in_format' not in kwargs.keys():
            in_format = self._extensions[filename.rsplit('.',1)[-1].lower()]
        else:
            in_format = kwargs['in_format']

        self.parser = self._parser_methods[in_format]()
        self.data = self.parser.parseFile(filename)
        '''try:
            self.parser = self._parser_methods[in_format]()
            self.data = self.parser.parseFile(filename)

        except:
            print("input file format not supported")'''

    def write(self, filename, out_format = 'Vamas'):
        self.writer = self._write_methods[out_format]()
        data = self.data
        self.writer.write(data, filename)
        '''try:
            self.writer = self._write_methods[out_format]()
            data = self.data
            self.writer.write(data, filename)
        except:
            print("output format not supported")'''
