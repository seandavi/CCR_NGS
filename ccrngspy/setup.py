from setuptools import setup, find_packages
import sys, os

version = '0.0'

setup(name='ccrngspy',
      version=version,
      description="'Python code to support ngs at NCI/CCR'",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='',
      author_email='seandavi@gmail.com',
      url='http://github.com/seandavi/CCR_NGS',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
