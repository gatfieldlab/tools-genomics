# gatlab-tools-genomics
This package provides a collection of data structures and utility functions for genomics-related scripts

## Installation

Using `pip` with our Python Package Server as an extra index:

    pip install gatlab-tools-genomics --extra-index-url https://gatfieldlab.github.io/python-package-server/
    
From `setup.py` with `setuptools`:

    install_requires = ['gatlab-tools-genomics-VERSION']
    dependency_links = ['https://github.com/gatfieldlab/tools_genomics/tarball/master#egg=gatlab-tools-genomics-VERSION']

Replace `VERSION` with the current version, for example 0.1.0.

To use a specific version from previous releases:

    dependency_links = ['https://github.com/gatfieldlab/tools_genomics/archive/v0.1.0.tar.gz#egg=gatlab-tools-genomics-0.1.0']
    
Support for `dependency_links` and `setuptools` will not be continued and will be removed in the next releases. Currently, such 
installs can only work with `python setup.py develop` from a local copy of the repo.

## Usage

Modules should be imported using the `gatlab_tools` namespace:

    from gatlab_tools.genomics import structure

