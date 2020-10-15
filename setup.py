#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("Pypi_README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

# with open("requirements.txt") as f:
#     requirements = f.read().splitlines()

# with open("requirements_dev.txt") as f:
#     requirements_dev = f.read().splitlines()

requirements = [
    "certifi==2020.6.20",
    "cycler==0.10.0",
    "fastcache==1.1.0",
    "kiwisolver==1.2.0",
    "matplotlib==3.3.2",
    "mpmath==1.1.0",
    "numpy==1.19.2",
    "pyparsing==2.4.7",
    "python-dateutil==2.8.1",
    "scipy==1.5.2",
    "six==1.15.0",
    "sympy==1.6.2",
    "tornado==6.0.4",
    "Click==7.1.2",
    "Sphinxcontrib-bibtex",
]
requirements_dev = [
    "pip==20.2.3",
    "bump2version==1.0.1",
    "wheel==0.35.1",
    "flake8==3.8.4",
    "tox==3.20.1",
    "coverage==5.3",
    "Sphinx==3.2.1",
    "twine==3.2.0",
]

setup_requirements = []

test_requirements = ["pytest>=3"]

setup(
    author="Zachary Streeter",
    author_email="zlstreeter@lbl.gov",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Exterior Complex Scaling Finite-Element Method Discrete Variable Representation grid for general physics problems.",
    entry_points={
        "console_scripts": [
            "quantumgrid=quantumgrid.cli:main",
            "ecs_femdvr_time_dep_h2=quantumgrid.ECS_FEMDVR_diatomic_time_dep_vibration_H2:main",
            "ecs_femdvr_time_indep_h2=quantumgrid.ECS_FEMDVR_diatomic_time_indep_vibration_H2:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords="quantumgrid",
    name="quantumgrid",
    packages=find_packages(
        include=["quantumgrid", "femdvr.py", "potential.py"]
    ),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/zstreeter/quantumgrid",
    version="0.1",
    zip_safe=False,
)
