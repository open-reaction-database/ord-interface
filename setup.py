#! /usr/bin/env python
# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Installer script."""
import setuptools

with open("README.md") as f:
    long_description = f.read()


if __name__ == "__main__":
    setuptools.setup(
        name="ord-interface",
        description="Interface for the Open Reaction Database",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/Open-Reaction-Database/ord-interface",
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: Apache Software License",
            "Operating System :: OS Independent",
        ],
        packages=setuptools.find_packages(),
        package_data={
            "ord_interface.visualization": [
                "template.html",
                "template.txt",
            ],
        },
        python_requires=">=3.10",
        install_requires=[
            "docopt>=0.6.2",
            "flask>=1.1.2",
            "flask-talisman>=1.0.0",
            "jinja2>=2.0.0",
            "ord-schema==0.3.71",
            "pandas>=1.0.4",
            "protobuf==4.22.3",
            "rdkit>=2021.9.5",
            "psycopg2>=2.8.5",
            "pygithub>=1.51",
            "requests>=2.24.0",
        ],
        extras_require={
            "tests": [
                "black>=22.3.0",
                "pylint>=2.13.9",
                "pytest>=6.2.5",
                "pytype>=2022.5.19",
            ],
        },
    )
