from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='taxamatcher',
      version='1.0',
      description='Easy to use taxonomy correlator',
      long_description=readme(),
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='taxonomy taxon taxa lineage taxonomy-assignment greengenes silva taxonomic-classification  bioinformatics metagenomics',
      url='https://github.com/mmtechslv/taxamatcher',
      author='Farid MUSA',
      author_email='farid.musa.h@gmail.com',
      license='Apache',
      packages= find_packages(),
      install_requires=['click','pandas','numpy'],
      scripts=['scripts/taxamatcher'],
      zip_safe=False)