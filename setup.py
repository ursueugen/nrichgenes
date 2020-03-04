import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nrich-ursueugen", # Replace with your own username
    version="0.0.1",
    author="Eugen Ursu",
    author_email="ursu_eugen@hotmail.com",
    description="Small package for gene enrichment, for internal use.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ursueugen/nrichgenes",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5'
)