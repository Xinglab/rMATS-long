import subprocess

import setuptools

import Cython.Build

COMMIT_FILE_NAME = 'source_code_commit.txt'


def get_source_code_commit():
    command = ['git', 'log', '--pretty=format:%H', '-n', '1']
    try:
        process = subprocess.run(command, capture_output=True, check=True)
    except:
        return 'unknown'

    commit_id = process.stdout.decode()
    return commit_id


def record_source_code_commit():
    commit = get_source_code_commit()
    commit_file_path = 'src/rmats_long/{}'.format(COMMIT_FILE_NAME)
    with open(commit_file_path, 'wt') as handle:
        handle.write('{}\n'.format(commit))


record_source_code_commit()

setuptools.setup(
    packages=['rmats_long'],
    package_dir={'': 'src'},
    package_data={'rmats_long': ['*.R', COMMIT_FILE_NAME]},
    ext_modules=Cython.Build.cythonize('src/rmats_long/*.py'),
)
