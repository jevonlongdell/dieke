import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dieke",
    version="0.0.1",
    author="Jevon Longdell",
    author_email="jevon.longdell@gmail.com",
    description="Crystal field calculations for the rare earths",
    long_description=long_description,
#   long_description_content_type="text/markdown",
    url="https://bitbucket.org/jevonlongdell/spectrapy",
    packages=['dieke'],
    install_requires=[
        'pandas',
        'xlrd',
        'numpy',
        'scipy',
        'matplotlib'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ),
    include_package_data=True
)
