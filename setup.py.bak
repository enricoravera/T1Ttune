#! /usr/bin/env python3

from setuptools import setup
from pathlib import Path

this_dir = Path(__file__).parent
long_description = (this_dir / "README.md").read_text(encoding="utf-8")

setup(
        name='t1ttune',
        version="0.0.1a",
        author='Enrico Ravera, Francesco Bruno, Letizia Fiorucci',
        author_email='ravera@cerm.unifi.it',
        description='A collection of scripts for setting up NMR protein dynamics experiments.',
        url='https://github.com/enricoravera/t1ttune',
        long_description=long_description,
        long_description_content_type='text/markdown',
        classifiers=[
            'Programming Language :: Python :: 3',
            'Operating System :: OS Independent',
            'License :: OSI Approved :: GPL-3.0 License',
            ],
        license='LICENSE.txt',
        install_requires=['numpy',
                          'scipy',
                          'lmfit',
                          'klassez',
                          'matplotlib>=3.8',
                          ],
        extras_require=None,
        packages=['t1ttune'],
        include_package_data=True,
        python_requires='>=3.10',
        )
