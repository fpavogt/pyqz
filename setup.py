from setuptools import setup, find_packages  # Always prefer setuptools over distutils

setup(
    name='pyqz',
    version='0.8.4.0',
    author='F.P.A. Vogt',
    author_email='frederic.vogt@alumni.anu.edu.au',
    packages=['pyqz',],
    url='http://fpavogt.github.io/pyqz/',
    download_url='https://github.com/fpavogt/pyqz/releases/tag/v0.8.4',
    license='GNU General Public License',
    description='Python link to MAPPINGS grids of HII regions.',
    long_description=open('README').read(),
        install_requires=[
        "numpy >= 1.12.1",
        "scipy >= 0.19.0",
        "matplotlib >= 2.0.1",
    ],
    include_package_data=True, # So that non .py files make it onto pypi, and then back !
)