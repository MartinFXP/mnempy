from setuptools import setup, find_packages

with open('README.md') as readme_file:
    README = readme_file.read()

with open('HISTORY.md') as history_file:
    HISTORY = history_file.read()

setup_args = dict(
    name='mnempy',
    version='0.9.0',
    description='Method to compute a mixture of nested effects models',
    long_description_content_type="text/markdown",
    long_description=README + '\n\n' + HISTORY,
    license='MIT',
    packages=find_packages(),
    author='Martin Pirkl',
    author_email='martinpirkl@yahoo.de',
    keywords=['Network', 'Perturbation', 'Bayesian'],
    url='https://github.com/MartinFXP/mnempy',
    download_url='https://pypi.org/project/mnempy/'
)

install_requires = [
    'numpy'
]

if __name__ == '__main__':
    setup(**setup_args, install_requires=install_requires)