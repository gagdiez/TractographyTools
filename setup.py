''' Setup file '''
from os import path
from distutils import extension

from setuptools import setup
import numpy

package_name = 'tractools'
cli_module = package_name + '.cli'
utils_module = package_name + '.utils'

setup(name=package_name,
      version='0.3',
      description='Tools to simplify the process of making tractography',
      url='http://github.com/gagdiez/tractools',
      author='Gallardo Diez, Guillermo Alejandro',
      author_email='gallardo@cbs.mpg.de',
      include_package_data=True,
      packages=[cli_module, utils_module],
      scripts=['scripts/tt_csd', 'scripts/tt_seeds_from_labeled_volume',
               'scripts/tt_tractography'],
      zip_safe=False)
