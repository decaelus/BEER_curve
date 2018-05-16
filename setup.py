from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='BEER_curve',
      version='0.1',
      author='Brian Jackson',
      author_email='bjackson@boisestate.edu',
      classifiers=[
        'Development Status :: 1 - Planning',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
      ],
      license=['GNU GPLv3'],
      packages=['BEER_curve'],
      install_requires=['PyAstronomy'],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      zip_safe=True)
