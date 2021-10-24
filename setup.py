from setuptools import setup, find_packages

setup(
    name = 'pytunnel', 
    version = '0.0.1', 
    description="Python implementation of the Crank-Nicholson Method for solving the Time-Dependent Schrodinger Equation for a Gaussian wave packet tunneling through a potential barrier.",
    url='https://gitlab.au.dk/ausa/oliskir/pytunnel',
    author='Oliver Kirsebom',
    author_email='oliver.kirsebom@gmail.com',
    license='GNU General Public License v3.0',
    install_requires=[
      'numpy',
      'scipy',
      'matplotlib',
      'Pint',
    ],
    packages=find_packages(),
    entry_points = {"console_scripts": ["pytunnel = pytunnel.pytunnel_script:main"]},
#    scripts=['pytunnel/bin/pytunnel_script.py'],
)
