import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

setuptools.setup(
    name="extendedMD",
    version="1",
    author="Maria Ines Silva",
    author_email="misilva73@gmail.com",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        'Documentation': 'https://extendedmd.readthedocs.io',
        'Source': 'https://github.com/misilva73/extendedMD',
    },
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    python_requires='>=3',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
