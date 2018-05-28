from setuptools import setup

# read the contents of your README file
# https://packaging.python.org/guides/making-a-pypi-friendly-readme/
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='BEER_curve',
    version='0.2',
    description='A very small package to model the BEaming, Ellipsoidal variations, and Reflected/emitted light from low-mass companions',
    author='Brian Jackson',
    author_email='bjackson@boisestate.edu',
    url='https://github.com/decaelus/BEER_curve',
    download_url = 'https://github.com/decaelus/BEER_curve/archive/0.1.tar.gz',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
    ],
    license=['GNU GPLv3'],
    packages=['BEER_curve'],
    install_requires=['PyAstronomy', 'statsmodels'],
    long_description=long_description,
    long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
    zip_safe=True)
