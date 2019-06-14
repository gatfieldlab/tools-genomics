from setuptools import setup

setup(
    name='gatlab-tools-genomics',
    version='0.1.0',
    description='Collection of data structures and '
                'utility functions for genomics',
#   url='NA',
    author='Bulak Arpat',
    author_email='bulak.arpat@gmail.com',
    license='GPLv3',
    packages=['gatlab_tools.genomics'],
    entry_points={
        'console_scripts': ['region = genomics.region:main']
    },
    install_requires=['numpy', 'gatlab-tools-extended==0.1.0'],
    dependency_links=['https://github.com/gatfieldlab/tools_extended#egg=gatlab-tools-extended-0.1.0'],
    zip_safe=False
)
