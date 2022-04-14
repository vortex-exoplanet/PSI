# Import all submodules.
from . import aperture
from . import psi

# Import all core submodules in default namespace.
from .aperture import *
from .psi import *

# Export default namespaces.
__all__ = []
__all__.extend(aperture.__all__)
__all__.extend(psi.__all__)

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
