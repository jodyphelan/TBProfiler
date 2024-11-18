from .text import *
from .reformat import *
from .collate import *
from .utils import *
from .spoligotyping import *
from .output import *
from .snp_dists import *
from .docx import *
from abc import ABC, abstractmethod

__version__ = "6.4.1"


class ProfilePlugin:
    """
    A class to define a plugin for tbprofiler
    """
    
    @abstractmethod
    def pre_process(self,args):
        """Generic pre-check method"""
        pass

    @abstractmethod
    def run(self):
        """Generic run method"""
        pass

    @abstractmethod
    def post_process(self,args):
        """Generic post-check method"""
        pass
