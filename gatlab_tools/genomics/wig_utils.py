#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provides tools to manipulate wig files
"""

import warnings
import numpy as np

ACCEPTED_MODES = ('w',)
DEF_TRACK_INFO = {'color': '255,0,0', 'maxHeightPixels': '80:60',
                  'graphType': 'bar', 'windowinfFunction': 'mean',
                  'coords': '1', 'scaleType': 'linear', 'type': 'wiggle_0',
                  'featureVisibilityWindow': '-1', 'gffTags': 'off',
                  'autoScale': 'on'}


class WigtoolsException(Exception):
    """
    A simple custom Exception class for wigtools.py
    """
    def __init__(self, m):
        message = "[wigtools] {}".format(m)
        super(WigtoolsException, self).__init__(message)


class WigFile(object):
    """
    A class for easy handling of wig files"
    """
    def __init__(self, file_path, mode, track_info=None, **kwargs):
        if not isinstance(file_path, str):
            raise WigtoolsException(
                'WigFile arg 1 has to be a string of file path')
        if mode not in ACCEPTED_MODES:
            raise Exception(
                'WigFile arg 2 has to one of {}'.format(ACCEPTED_MODES))
        # Private inits
        self.__file_path = file_path
        self.__mode = mode
        self.__track_info = None
        self.__handle = None
        self.__info_written = False
        self.__track_line = None
        self.__cur_step_type = None
        self.__last_pos = None
        # Public inits
        self.track_info = track_info
        super().__init__(**kwargs)

    def __del__(self):
        self.__handle.flush()
        self.__handle.close()

    # Private methods

    def __open_file(self):
        self.__handle = open(self.__file_path, self.__mode)

    def __writeln(self, line):
        if not self.is_open:
            raise WigtoolsException('WigFile object file is not opened')
        if self.mode not in ['w']:
            raise WigtoolsException('WigFile object is not in a writable mode')
        self.__handle.write("{}\n".format(line))

    # Attribute properties

    @property
    def track_info_written(self):
        """Returns if the track info line was already written or not"""
        return self.__info_written

    @property
    def mode(self):
        """Returns the file mode of the WigFile object"""
        return self.__mode

    @property
    def track_info(self):
        """Returns the track information of the WigFile object"""
        return self.__track_info

    @track_info.setter
    def track_info(self, value):
        if value is None:
            return
        if isinstance(value, dict):
            self.__track_info = value
        else:
            raise Exception(
                'WigFile track_info arg, if set, has to a dict of track info')

    @property
    def is_open(self):
        """
        Returns if the WigFile object file is opened or not
        """
        return self.__handle is not None

    @property
    def track_info_ok(self):
        """
        Returns if the track info meets minimal requirements or not
        """
        try:
            result = self.__track_info['type'] == 'wiggle_0'
        except (TypeError, KeyError):
            result = False
        return result

    # Public methods

    def write_trackinfo(self, track_info=None):
        """
        Writes the track info as a first line to the wig file
        """
        if self.mode == 'w' and not self.track_info_written:
            if not self.is_open:
                self.__open_file()
            if not self.track_info:
                self.track_info = track_info or DEF_TRACK_INFO
            if not self.track_info_ok:
                raise WigtoolsException(
                    'Minimal track info requirements are not met')
            track_line = 'track ' + ' '.join(['{}={}'.format(k, v) for k, v in
                                              self.track_info.items()])
            self.__writeln(track_line)
            self.__info_written = True
            self.__track_line = track_line
        elif self.mode == 'w' and self.track_info_written:
            warnings.warn(
                'Track info line has already been written!', RuntimeWarning)
        else:
            warnings.warn(
                "Track info line can only be written in 'w' mode",
                RuntimeWarning)

    def write_new_chrom(self, name, step_type,
                        start=None, step=None, span=None):
        """
        Writes a new chromosome line with the given settings. Also uses the
        same settings for writing later any data lines that will be passed on
        """
        if not (isinstance(step_type, str) or
                step_type in ('fixed', 'variable')):
            raise WigtoolsException(
                "step_type argument must be either 'fixed' or 'variable'")
        if not (isinstance(name, str) or name == ''):
            raise WigtoolsException(
                'name argument has be a non-empty string')
        chrom_line = 'chrom={}'.format(name)
        if step_type == 'fixed':
            if not (isinstance(start, int) and isinstance(step, int)):
                raise WigtoolsException(
                    'for fixedStep, start and step arguments have to be int')
            if start < 1 or step < 1:
                raise WigtoolsException(
                    'start and step arguments have to be positive')
            chrom_line = "fixedStep {} start={} step={}".format(
                chrom_line, start, step)
            if span is not None:
                if not isinstance(span, int) or span < 1:
                    raise WigtoolsException(
                        'span, is set, has to be a positive integer')
                chrom_line += ' span={}'.format(span)
        if step_type == 'variable':
            chrom_line = "variableStep {}".format(chrom_line)
            if span is not None:
                if not isinstance(span, int) or span < 1:
                    raise WigtoolsException(
                        'span, is set, has to be a positive integer')
                chrom_line += " span={}".format(span)
        if not self.track_info_written:
            warnings.warn(
                """Can not write a new chromosome when track info has not been
                written yet. Enforcing the writing of the track info line...
                """)
            self.write_trackinfo()
        self.__writeln(chrom_line)
        self.__cur_step_type = step_type
        self.__last_pos = 0

    def write_data(self, data, start=None):
        """
        Writes the data conforming the current step type. This means,
         - if step type is variable:
            - if data is 1D array, it will be interpreted as a continous array
                starting from start argument if set, otherwise from last pos
            - if data is 2D array, it will be written as it is, ignoring any
                value start may have
        - if step type is fixed:
            data will be interpreted as a continuous array, start will be
            ignored
        """
        last_pos = start or self.__last_pos
        if isinstance(data, (int, float)):
            data = ([last_pos+1], [data])
        else:
            if isinstance(data, (list, tuple)):
                np_data = np.array(data)
            else:
                np_data = data
            if not isinstance(np_data, np.ndarray):
                raise WigtoolsException(
                    ('data argument has to be one of: int, ' +
                     'float, tuple, list or a numpy array'))
            if not issubclass(np_data.dtype.type, (np.integer, np.floating)):
                raise WigtoolsException(
                    'data has to be an array-like of integers or floats')
            if len(np_data.shape) == 1:
                data = (range(last_pos+1, last_pos+np_data.shape[0]+1),
                        np_data.tolist())
            else:
                if len(np_data.shape) > 2 or np_data.shape[0] != 2:
                    raise WigtoolsException(
                        'data can be max a 2D array with 2 rows')
                if self.__cur_step_type == 'fixed':
                    raise WigtoolsException(
                        'a 2D array can only be written for variableStep')
                if not np.all(np.mod(np_data[0], 1) == 0):
                    raise WigtoolsException(
                        ('first row of 2D data structure must contain ' +
                         'integers or whole numbers only - for positions'))
                data = (np_data[0].astype(int).tolist(), np_data[1].tolist())
            if self.__cur_step_type == 'fixed':
                for value in data[1]:
                    self.__writeln(value)
            elif self.__cur_step_type == 'variable':
                for i in range(len(data[0])):
                    self.__writeln("{} {}".format(data[0][i], data[1][i]))
            self.__last_pos = data[0][-1]
            print(self.__last_pos)
