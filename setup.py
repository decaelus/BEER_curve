from setuptools import setup

def readme():
    with open('README.rst') as f:
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
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      zip_safe=True)
