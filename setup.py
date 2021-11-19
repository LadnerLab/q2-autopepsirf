from setuptools import setup, find_packages


setup(
    name="q2-autopepsirf",
    version='0.0.1.dev',
    packages=find_packages(),
    package_data={},
    author="Annabelle Brown",
    author_email="annabelle811@live.com",
    description="Auto-Run q2-pepsirf and q2-ps-plot",
    license='Apache-2.0',
    url="https://github.com/LadnerLab/q2-autopepsirf",
    entry_points={
        'qiime2.plugins': ['q2-autopepsirf=q2_autopepsirf.plugin_setup:plugin']
    },
    zip_safe=False,
)