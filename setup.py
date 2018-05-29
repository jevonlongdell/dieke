import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="spectrapy",
    version="0.0.1",
    author="Jevon Longdell",
    author_email="jevon.longdell@gmail.com",
    description="Crystal field calculations for the rare earths",
    long_description=long_description,
#   long_description_content_type="text/markdown",
    url="https://bitbucket.org/jevonlongdell/spectrapy",
    packages='spectrapy',
    install_requires=['numpy','pandas'],
    zip_safe=False)
