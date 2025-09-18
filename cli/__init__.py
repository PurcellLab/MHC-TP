__version__ = "1.1.1-beta"


from warnings import filterwarnings
from rich.traceback import install

install(show_locals=True)  # type: ignore

from .HLAfreq import *