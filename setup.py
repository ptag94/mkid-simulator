# setup.py
from setuptools import setup, find_packages

setup(
    name="mkid-simulator",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        # Liste des dépendances, par exemple : 'requests', 'numpy'
    ],
    author="TAGNON Paul",
    author_email="ptag94@gmail.com",
    description="Package for Micro-wave Kinetic Inductance Detector simulation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ptag94/mkid-simulator",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
