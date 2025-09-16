import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(

    #Basic description
    name='TSpy',
    author="Deivy Mercado",
    author_email="david231097@gmail.com",
    description="Tidal signal with Python and SPICE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DeivyMercado/TSpy",
    keywords='astrodynamics geophysics tides spiceypy',
    license='MIT',

    #Clasifier
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
        ],
    version='0.0.1',

    #Files
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    
    # ######################################################################
    # ENTRY POINTS
    # ######################################################################
    #entry_points={
        #'console_scripts': ['install=pymcel.install:main'],
    #},

    # ######################################################################
    # TESTS
    # ######################################################################
    #test_suite='nose.collector',
    #tests_require=['nose'],

    #Dependencies
    install_requires=['numpy',
                      'scipy',
                      'spiceypy'],

    #Options
    include_package_data=True,
    package_data={'tspy': ['data/spice_kernels.json']
                  },
)
