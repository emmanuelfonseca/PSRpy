PSRpy
=====

This repo is an evolving collection of Python3 functions and modules for pulsar and radio-transient science. This is a work in progress, and documentation will slowly be added over time. However, the basic aspects of `PSRpy` are described below.

## Installation

`PSRpy` can be installed by cloning the repo and running `pip` in the following way:

``` 
pc> git clone https://github.com/emmanuelfonseca/PSRpy.git
pc> cd PSRpy
pc> pip install .
```

## Usage

Once installed, `PSRpy` can be imported into a Python3 interpreter and its utilities accessed in the standard package format. For example, if you wish to compute the orbital elements of a binary system based on observed parameters, you can obtain information on the avaialble, relevant functions in the following manner:

``` python
pc> python
>>> import PSRpy.orbit.elements as elem
>>> help(elem)
```

There are a growing number of Python3 scripts that use the `PSRpy` package in different ways and context. These scripts are located in the `bin/` subdirectory of the `PSRpy` package. Feel free to suggest or submit modifications and/or contribute scripts you think would be of general interest!
