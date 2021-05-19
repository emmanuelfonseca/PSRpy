One of the common operations in pulsar astronomy is the reading, writing, and editing of timing solutions or "parameter files" (colloquially referred to as "parfiles"). The `PSRpy` package contains a `Parfile` object that centralizes operations related to timing solutions and their manipulation.

## Instantiating the Parfile() Object
The `Parfile` object currently does not depend on a user-supplied file, and can be instantiated in the following way:

``` python
> python
>>> from PSRpy.parfile import Parfile
>>> input_parfile = Parfile()
>>>
```

## Format of Parfile Attributes
The `Parfile` object mostly consists of attributes that correspond to timing-model parameters used in community-developed timing software, such as `tempo`, `tempo2`, and `PINT`. With the exception of TOA-noise terms, all `Parfile` attributes are Python dictionaries with the following format:

``` python
>>> pmra = getattr(input_parfile, "PMRA")
>>> print(pmra)
{'value': None, 'flag': None, 'error': None}
```

Nearly all fit parameters -- ones subject to least-squares fitting -- are listed in parfiles in "NAME VALUE FLAG ERROR" format; the `Parfile` object organizes these data for all parameters -- whether they are fit or configuration parameters -- into dictionaries to make subsequent access and manipulation Pythonic.

## Reading Models from a File
The `Parfile` object also contains several methods for I/O and editing of timing solutions in ASCII format. To read a parfile and overload all relevant attributes in the above `Parfile` object, use the `read()` method:

``` python
>>> input_parfile.read("B1534+12.par")
>>> pmra = getattr(input_parfile, "PMRA")
>>> print(pmra)
{'value': 1.482, 'flag': 0, 'error': 0.007}
```

In the above example, the "B1534+12.par" ASCII timing model is taken from the public ATNF pulsar catalog.

## Editing Attribute Data
Depending on the context, users may need to update fit flags, change parameter values, or even introduce new parameters not supplied in the input timing solution. The `set()` method in the `Parfile` object allows you to make such changes by supplying a dictionary containing data to be loaded into the object.

For example, let's assume we want to enable fitting of the `PMRA` parameter; this opertation means we need to change the value of the `flag` key for the `Parfile` attribute from `0` to `1`. Using the `set()` method, we need to supply a dictionary with one key/value pair:

``` python
>>> print(pmra)
{'value': 1.482, 'flag': 0, 'error': 0.007}
>>> input_parfile.set("PMRA", {"flag": 1})
>>> pmra = getattr(input_parfile, "PMRA")
>>> print(pmra)
{'value': 1.482, 'flag': 1, 'error': 0.007}
```

Notice that only the value of the `flag` key changes with the above `set()` call. If we instead wish to set the `PMRA` value fixed to zero -- effectively neglecting proper motion along the right ascension coordinate -- we merely need to supply the right dictionary:

``` python
>>> print(pmra)
{'value': 1.482, 'flag': 1, 'error': 0.007}
>>> input_parfile.set("PMRA", {"value": 0., "flag": 0, "error": None})
>>> pmra = getattr(input_parfile, "PMRA")
>>> print(pmra)
{'value': 0.0, 'flag': 0, 'error': None}
```

## Writing Parfiles
The `write()` method will examine all `Parfile` attributes and write those with initialized values to disk:

``` python
>>> input_parfile.write(outfile="new.par")
```

A new timing solution should appear in your local working area as the "new.par" ASCII file.

## Creating Models from Scratch
As discussed above, the `Parfile` object can work independently of external information. An example is the instance where a user would like to simulate timing models and their corresponding TOAs based on custom simulation software. The following script demonstrates how the `Parfile` object and its methods can be used to create a fake parfile.

``` python
#! /usr/bin/env python
  
from PSRpy.parfile import Parfile

# initialize Parfile object.
fake_par = Parfile()

# load in fake data.
fake_par.set("PSR", {"value": "J1234+5678"})
fake_par.set("RAJ", {"value": "12:34:00.0"})
fake_par.set("DECJ", {"value": "56:78:00.0"})
fake_par.set("F0", {"value": 0.3333333333})
fake_par.set("F1", {"value": -1.0e-15})
fake_par.set("PEPOCH", {"value":53005.0})

# write to file.
fake_par.write(outfile="fake.par")
```
