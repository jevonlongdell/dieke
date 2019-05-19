import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dieke",
    version='0.3.0',
    author="Jevon Longdell",
    author_email="jevon.longdell@gmail.com",
    description="Crystal field calculations for the rare earths",
    long_description=long_description,
    url="https://github.com/jevonlongdell/dieke",
    packages=['dieke'],
    install_requires=[
        'pandas',
        'xlrd',
        'numpy',
        'scipy',
        'matplotlib',
        'sympy'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2",
        "Operating System :: OS Independent",
    ),
    include_package_data=True
)
