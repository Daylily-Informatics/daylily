from setuptools import setup, find_packages

setup(
    name="daylily",
    version="0.7.196t",
    packages=find_packages(),
    install_requires=[
        # Add dependencies here
    ],
    package_data={
        "": ["scripts/*.sh"],  # Include all `.sh` files under the `scripts/` directory
    },
    entry_points={
        "console_scripts": [
            "calc_daylily_aws_cost_estimates=bin.calc_daylily_aws_cost_estimates:main",  # Maps "calc_costs" command to Python script
        ],
    },
)
