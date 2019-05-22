#!/usr/bin/env python

import logging
import os


def fetch_abspath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def set_log_config(debug=None, info=None):
    if debug:
        lv = logging.DEBUG
    elif info:
        lv = logging.INFO
    else:
        lv = logging.WARNING
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=lv
    )
