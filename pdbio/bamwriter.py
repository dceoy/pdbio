#!/usr/bin/env python
"""
BAM File Writer.
https://github.com/dceoy/pdbio
"""

import logging
import os
import subprocess


class BamFileWriter(object):
    def __init__(self, out_bam_path, samtools, n_thread=1):
        self.__logger = logging.getLogger(__name__)
        self.__bam_path = out_bam_path
        self.__logger.debug('Write STDIN into BAM: {}'.format(self.__bam_path))
        self.__samtools = samtools
        args = [
            self.__samtools, 'view', '-@', str(n_thread), '-bS', '-', '-o',
            out_bam_path
        ]
        self.__logger.debug('STDIN => `{}`'.format(' '.join(args)))
        self.proc = subprocess.Popen(
            args=args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

    def write(self, *args, **kwargs):
        self.proc.stdin.write(*args, **kwargs)
        self.proc.stdin.flush()

    def close(self):
        self.__logger.debug('Finish writing BAM: {}'.format(self.__bam_path))
        self.proc.stdin.close()
        self.proc.wait()
        if self.proc.returncode == 0:
            subprocess.run([self.__samtools, 'quickcheck', self.__bam_path])
        else:
            self.__logger.error(
                'STDERR from subprocess `{0}`:{1}{2}'.format(
                    self.proc.args, os.linesep,
                    self.proc.stderr.decode('utf-8')
                )
            )
            raise subprocess.CalledProcessError(
                returncode=self.proc.returncode, cmd=self.proc.args,
                output=self.proc.stdout, stderr=self.proc.stderr
            )

    def __enter__(self):
        return self

    def __exit__(self, ex_type, ex_value, trace):
        args = {'ex_type': ex_type, 'ex_value': ex_value, 'trace': trace}
        self.__logger.debug('with-statement exit: {}'.format(args))
        self.close()

    def __del__(self):
        self.close()
