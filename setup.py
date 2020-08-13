from setuptools import setup

setup(
    name='gatlab-tools-genomics',
    version='0.1.1',
    description='Collection of data structures and '
                'utility functions for genomics',
#   url='NA',
    author='Bulak Arpat',
    author_email='bulak.arpat@gmail.com',
    license='GPLv3',
    packages=['gatlab_tools.genomics'],
    entry_points={
        'console_scripts': ['region = gatlab_tools.genomics.region:main']
    },
    install_requires=['numpy',
                      'gatlab-tools-extended==0.1.0',
                      'gatlab-tools-specialized==0.1.0'],
    # dependency_links is not used by pip, only by setuptools if used as "python setup.py install"
    # it will be removed in the following releases
    dependency_links=['https://github.com/gatfieldlab/tools_extended/tarball/master#egg=gatlab-tools-extended-0.1.0',
                      'https://github.com/gatfieldlab/tools_specialized/tarball/master#egg=gatlab-tools-specialized-0.1.0'],
    zip_safe=False
)
