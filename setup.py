##################################################################
#                                                                #
# tSPICE: Tidal Signal with Python and SPICE                     #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)       #
##################################################################
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='tspice',
    author="Deivy Mercado, Jorge I. Zuluaga, Gloria Moncayo",
    author_email="david231097@gmail.com, jorge.zuluaga@udea.edu.co, gloria.moncayo@udea.edu.co",
    description="Tidal signal with Python and SPICE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DeivyMercado/TSPICE",
    keywords='astrodynamics geophysics tides spice',
    license='AGPL-3.0-only',

    # ######################################################################
    # CLASSIFIER
    # ######################################################################
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        #"License :: OSI Approved :: GNU Affero General Public License v3 (AGPLv3)",
        "Operating System :: OS Independent",
    ],
    version='0.0.2',

    # ######################################################################
    # FILES
    # ######################################################################
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    
    # ######################################################################
    # TESTS
    # ######################################################################
    test_suite='tests',
    tests_require=['pytest'],

    # ######################################################################
    # DEPENDENCIES
    # ######################################################################
    install_requires=[
        'numpy',
        'scipy',
        'spiceypy',
        'matplotlib',
        'pandas',
        'openpyxl'
    ],
    
    python_requires='>=3.7',

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={'tspice': ['data/*.*']},
)