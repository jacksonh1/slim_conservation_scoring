import os

from setuptools import find_packages, setup

setup(
    name="local_conservation_pipeline",
    # mandatory
    version="0.1",
    # mandatory
    author="Jackson Halpin",
    packages=find_packages(
        where='src',
    ),
    include_package_data=True,
    package_dir={"": "src"}
)
