from setuptools import setup, find_packages

# run setup for PSRpy.
setup(
    name="PSRpy",
    version="0.0",
    packages=find_packages(),
    scripts=[
        "bin/compute_Shapiro_grid.py",
        "bin/compute_position_grid.py",
        "bin/compute_orbital_phase_epochs.py",
        "bin/fit_orbit.py",
        "bin/plot_residuals.py",
    ],

    # metadata to display on PyPI.
    author="Emmanuel Fonseca",
    author_email="efonseca@physics.mcgill.ca",
    description="PSRpy: a collection of utilities related to pulsar astronomy.",
    keywords="pulsar astronomy orbits tempo timing"
)
