#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Rohan Weeden; Brett Buzzanga
# Copyright 2020, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Global logging configuration
Inspired By:
https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
"""
import logging
import os, os.path as op
from logging import FileHandler, Formatter, StreamHandler


def log2file(fname='error.log'):
    """ Also write log to file e.g.: bbLogger.log2file(PATH) """
    fname  = op.join(fname, 'error.log') if op.isdir(fname) else fname
    errorfile_handler = FileHandler(fname)
    errorfile_handler.setFormatter(Formatter(
        "[{asctime}] {funcName:>20}:{lineno:<5} {levelname:<10} {message}",
        style="{", datefmt='%Y-%m-%d %H:%M:%S'))
    try:
        errorfile_handler.setLevel(logger.getEffectiveLevel())
    except:
        errorfile_handler.setLevel(logging.INFO)

    ## addded recently might not work
    # Check if the logger already has handlers to avoid duplication
    if not logger.hasHandlers():
        logger.addHandler(errorfile_handler)
        
    # To avoid logger propagation
    logger.propagate = False 


class UnixColorFormatter(Formatter):
    ## seems to work for piping output with gattaca
    blue     = '\033[94m'
    yellow   = "\x1b[33;21m"
    red      = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset    = "\x1b[0m"

    COLORS = {
        logging.INFO    : blue,
        logging.WARNING : yellow,
        logging.ERROR   : red,
        logging.CRITICAL: bold_red
    }


    def __init__(self, fmt=None, datefmt=None, style="%", use_color=False):
        super().__init__(fmt, datefmt, style)
        # Save the old function so we can call it later
        self.__formatMessage = self.formatMessage
        if use_color:
            self.formatMessage = self.formatMessageColor


    def formatMessageColor(self, record):
        message = self.__formatMessage(record)
        color = self.COLORS.get(record.levelno)

        if color:
            # message = "".join([color, *[i for i in message], self.reset])
            message = "".join([color,  message, self.reset])
        return message


class CustomFormatter(UnixColorFormatter):
    """Adds levelname prefixes to the message on warning or above."""

    def formatMessage(self, record):
        message = super().formatMessage(record)
        if record.levelno >= logging.DEBUG:
            message = ": ".join((record.levelname, message))
        return message


logger = logging.getLogger('cariSLR_git')
logger.setLevel(logging.INFO)

stdout_handler = StreamHandler(os.sys.stdout)
fmt     = '  {funcName}:{lineno:<6} [{asctime}] {message}'
datefmt = '%b %d, %H:%M:%S'

# Obj = CustomFormatter(fmt, datefmt, '{', use_color=os.name != "nt")
LogObj = CustomFormatter(fmt, datefmt, '{', use_color=os.sys.stdout.isatty())
stdout_handler.setFormatter(LogObj)

# stdout_handler.setFormatter(Formatter(
#     '{funcName:>10}:{lineno:<6} [{asctime}] {levelname:<10} {message}', style='{', datefmt='%b %d, %H:%M:%S'))

logger.addHandler(stdout_handler)
