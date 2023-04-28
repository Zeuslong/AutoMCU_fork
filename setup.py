from importlib.metadata import entry_points
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="amcu",
    version="0.0.1",
    description=("A package to quantify endmember fractions"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["rasterio", "pydantic", "pytest", "spectral", "numpy", "scipy"],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    entry_points={"console_scripts": ["amcu = amcu.cli:main"]},
)
