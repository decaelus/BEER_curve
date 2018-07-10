from setuptools import setup
import BEER_curve

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='BEER_curve',
    version=BEER_curve.__version__,
    description='A very small package to model the BEaming, Ellipsoidal variations, and Reflected/emitted light from low-mass companions',
    author='Brian Jackson',
    author_email='bjackson@boisestate.edu',
    url='https://github.com/decaelus/BEER_curve',
    download_url =
    'https://github.com/decaelus/BEER_curve/archive/'+BEER_curve.__version__+'.tar.gz',
    long_description = 'A very small package to model the BEaming, Ellipsoidal variations, and Reflected/emitted light from low-mass companions',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
    ],
    license=['MIT'],
    packages=['BEER_curve'],
    install_requires=['PyAstronomy', 'statsmodels'])
