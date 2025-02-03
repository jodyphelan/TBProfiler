import argparse
from pathogenprofiler.models import Variant
from .models import ProfileResult
from typing import List

class ProfilePlugin:
    """
    A class to define a plugin for tbprofiler
    """
    
    def pre_process(self,args):
        """Generic pre-check method"""
        pass

    def run(self):
        """Generic run method"""
        pass

    def post_process(self,args):
        """Generic post-check method"""
        pass

    def process_variants(self,args: argparse.Namespace, variants: List[Variant]):
        """Generic variant processing method"""
        pass

    def process_result(self,args: argparse.Namespace, result: ProfileResult):
        """Generic result procesing method"""
        pass