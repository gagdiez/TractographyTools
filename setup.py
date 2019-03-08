''' Setup file '''
from os import path
from distutils import extension

from setuptools import setup
import numpy

package_name = 'tractools'
cli_module = package_name + '.cli'

setup(name=package_name,
      version='0.2',
      description='Tools to simplify the process of making tractography',
      url='http://github.com/gagdiez/tractools',
      author='Gallardo Diez, Guillermo Alejandro',
      author_email='gallardo@cbs.mpg.de',
      include_package_data=True,
      packages=[cli_module],
      scripts=['scripts/tt_csd'],
      zip_safe=False)
