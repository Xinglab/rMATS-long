import argparse
import importlib.resources
import sys

from rmats_long import rmats_long_utils


# argparse is used for the help message and for checking the script choices.
# Only the 1st argument is given to parse_args() to
# avoid errors from any "--" args being passed to other scripts.
# The arguments after script are taken directly from sys.argv.
def parse_args(script_names):
    parser = argparse.ArgumentParser(
        description='Run a script from the rmats_long package')
    parser.add_argument(
        'script',
        choices=script_names,
        help='The script to run. One of: {}'.format(script_names))
    parser.add_argument('arguments',
                        nargs='*',
                        help='Arguments to pass along to the script')

    args_to_parse = sys.argv[1:2]
    remaining = sys.argv[2:]
    args = parser.parse_args(args_to_parse)
    args.arguments = remaining
    return args


def get_script_names():
    names = list()
    module_dir = importlib.resources.files('rmats_long')
    for file_object in module_dir.iterdir():
        if not file_object.is_file():
            continue

        name = file_object.name
        is_py_name = name.endswith('.py') and name != '__init__.py'
        is_r_name = name.endswith('.R')
        if is_py_name or is_r_name:
            names.append(name)

    names.sort()
    return names


def main():
    script_names = get_script_names()
    args = parse_args(script_names)
    rmats_long_utils.run_package_script(args.script, args.arguments)


if __name__ == '__main__':
    main()
