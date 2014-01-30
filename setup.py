import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
__version = 'a0.0.0'
setup(
    name = "nampy",
    version = __version,
    packages = find_packages(),
    zip_safe = False,
    #scripts = [''],
    #setup_requires = [],
    #install_requires = ['myitem>=0.1'],
    extras_require = {},

    package_data = {'': ['*.txt', '*.html','LICENSE','README','data/*','examples/*py','rfiles/*']},

    author = "Brian J Schmidt",
    author_email = "brianjschmidt@gmail.com",
    description = "",
    license = "GPL V3.0",
    keywords = "systems biology networks omics",
    url = "",
    test_suite = "",
    long_description = "",
    )
    
