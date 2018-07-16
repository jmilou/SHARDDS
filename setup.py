from setuptools import setup, find_packages

setup(
    name='library',
    version='0.1',
    description='Tool to create a library based on the SHARDDS program',
    url='https://github.com/jmilou/SHARDDS',
    author='Julien Milli',
    author_email='jmilli@eso.org',
    license='MIT',
    keywords='reference star differential imaging',
    packages=find_packages(),
    install_requires=[
        'numpy', 'scipy', 'astropy', 'matplotlib'
    ],
    zip_safe=False
)
