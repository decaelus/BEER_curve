from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='BEER_curve',
      version='0.2',
      author='Brian Jackson',
      author_email='bjackson@boisestate.edu',
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
