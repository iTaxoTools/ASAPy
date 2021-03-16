from multiprocessing import Process

import tempfile
import shutil
import pathlib
from datetime import datetime

from . import asapc
from . import param
from . import params

class PartitionAnalysis():
    """
    Container for input/output of ASAP core.
    """

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state

    def __init__(self, file):
        """
        """
        self.file = file
        self.useLogfile = False
        self.target = None
        self.results = None
        # self.time_format = '%x - %I:%M%p'
        self.time_format = '%FT%T'
        self.param = param.ParamList(params.params)

    def fetch(self, destination):
        """
        Copy results as a new directory.
        """
        if self.results is None:
            raise RuntimeError('No results to fetch.')
        shutil.copytree(self.results, destination)

    def run(self):
        """
        Run the ASAP core with given params,
        save results to a temporary directory.
        """
        kwargs = self.param.as_dictionary()
        kwargs['file'] = self.file
        kwargs['logfile'] = self.useLogfile
        kwargs['time'] = datetime.now().strftime(self.time_format)
        if self.target is not None:
            kwargs['out'] = self.target
        print(kwargs)
        asapc.main(kwargs)
        self.results = self.target

    def launch(self):
        """
        Should always use a seperate process to launch the ASAP core,
        since it uses exit(1) and doesn't always free allocated memory.
        Save results on a temporary directory, use fetch() to retrieve them.
        """
        # When the last reference of TemporaryDirectory is gone,
        # the directory is automatically cleaned up, so keep it here.
        self._temp = tempfile.TemporaryDirectory(prefix='asap_')
        self.target = pathlib.Path(self._temp.name).as_posix()
        p = Process(target=self.run)
        p.start()
        p.join()
        if p.exitcode != 0:
            raise RuntimeError('ASAP internal error, please check logs.')
        # Success, update analysis object for parent process
        self.results = self.target
