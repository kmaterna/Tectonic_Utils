import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Tectonic_Utils", 
    version="0.0.8",
    author="Kathryn Materna",
    author_email="kmaterna@berkeley.edu",
    description="A small package of useful geophysics utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kmaterna/Tectonic_Utils",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
