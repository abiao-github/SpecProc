from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="SpecProc",
    version="1.0.0",
    author="Niu, Hubiao",
    author_email="your.email@example.com",
    description="A spectral processing tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/specproc",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        "matplotlib",
        "PyQt5",
        "configparser",
    ],
    entry_points={
        'console_scripts': ['specproc=src.main:main',],
    },
)