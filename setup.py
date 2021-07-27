from setuptools import setup, find_packages

setup(name='metaMIC',
      description='metaMIC: Reference-free Misassembly Identification and Correction of metagenomic assemblies',
      long_description_content_type = 'text/markdown',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
      ],
      url='https://github.com/Seny-l/metaMIC',
      author='Senying Lai',
      license='MIT',
      packages = find_packages(exclude=['docs', 'tests', 'examples']),
      package_data={'metaMIC' : ['*.py','*.sh']},
      entry_points={
            'console_scripts' : ['metaMIC=metaMIC.__main__:main'],
      },
      zip_safe=False
)

