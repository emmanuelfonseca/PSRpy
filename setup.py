from setuptools import setup, find_packages
import numpy as np

# run setup for PSRpy.
setup(
    name="PSRpy",
    version="0.0",
    include_dirs=[np.get_include()],
    packages=find_packages(),
    scripts=[
        "PSRpy/bin/compute_chisq_grid.py",
        "PSRpy/bin/compute_Shapiro_pdfs.py",
        "PSRpy/bin/compute_position_grid.py",
        "PSRpy/bin/compute_orbital_phase_epochs.py",
        "PSRpy/bin/fit_orbit.py",
        "PSRpy/bin/plot_residuals.py",
    ],

    # metadata to display on PyPI.
    author="Emmanuel Fonseca",
    author_email="efonseca@physics.mcgill.ca",
    description="PSRpy: a collection of utilities related to pulsar astronomy.",
    keywords="pulsar astronomy orbits tempo timing"
)
