from distutils.core import setup
from setuptools import find_packages
setup(
    name='medusa-cobra',
    description='Ensemble modeling using constraint-based reconstruction and analysis (COBRA)',
    long_description='Medusa extends various algorithms and analyses used in the COBRA field to ensemble-scale, and provides a high-level user interface for analysis, visualization, and interpreting simulation results.',
    version='0.3.0',
    author='Gregory Medlock',
    author_email='gmedlo@gmail.com',
    url='https://github.com/opencobra/medusa',
    download_url='https://github.com/opencobra/medusa/archive/v0.3.0.tar.gz',
    keywords=['computational biology','systems biology','biology'],
    license='MIT',
    classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    ],
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=['cobra>=0.13.0'],
    package_data={'':  ['test/data/*']}
)
